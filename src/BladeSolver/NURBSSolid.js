"use strict";
/******
Copyright 2000-2015 Mark Hampsey // escapee.org

Portions of this code are adapted from routines in 'The NURBS Book' by Piegl and Tiller:
http://www.springer.com/gp/book/9783642973857
*******/

class NURBSSolid{
	constructor(parameters){
		this.curves;		//2D array of SurfaceCurve objects, one for each direction
		this.q;   		//degree of surface in V direction
		this.r;   		//degree of surface in W direction	
		this.V;			//the knot vector in the V direction
		this.W;			//the knot vector in the W direction
		
		this.scale = 1.0;
		this.previousScale = 1.0;
		
		this.alignedAtPercentChord = undefined;
		this.centroidAtPercentChord = undefined;

		this.wDegree; //degree of solid in W direction
		this.vDegree; //degree of solid in V direction

		if(parameters)this.initialise(parameters);
	}

	//optional constructor
	initialise(parameters){
		//input: parameters{surfaceCurves: an array of SurfaceCurve objects,
		//					vDegree: degree of surface in V direction (int), 
		//					vBspline: a pre-defined Bspline curve representation of V geometry. (optional)
		//					wDegree, degree of surface in W direction (int), 
		//					alignedAtPercentChord: determines how spanwise (V) curves are aligned [0.,1.], 
		//					centroidAtPercentChord: determines how thickness curves (W) are aligned [0.,1.]}
		//surfaceCurves defines the backbone of the surface S(u,v)
	
		if(!parameters || !parameters.surfaceCurves){
			console.log("NURBSSolid: not enough info for model initialisation.");
			return;
		}	
		
		/////
		//V direction
		/////
		this.curves = []; //an empty array that will become an array of surfaceCurve arrays;
		this.curves.push(parameters.surfaceCurves);
		
		//default alignment of each curve to the leading edge.
		if(parameters.alignedAtPercentChord !== undefined && parameters.alignedAtPercentChord>=0.)
			this.alignedAtPercentChord = parameters.alignedAtPercentChord;
		
		//Set the origin of each curve to the chosen alignment point.
		for(var i=0;i<this.curves[0].length;i++)this.curves[0][i].moveOriginToPercentChord(this.alignedAtPercentChord);

		if(parameters.vBspline){
			//if V geometry has been precomputed as a Bspline curve, use that curve's degree and knot vector.
			this.vDegree = parameters.vBspline.p;
			this.q = this.vDegree;
			this.V = [];
			for(var i=0;i<parameters.vBspline.U.length;i++)this.V.push(parameters.vBspline.U[i]);
		}
		else{
			if(parameters.vDegree != undefined)this.vDegree = parameters.vDegree;		
			this.calculate_knots(0); //0 == create initial linear surface/solid in V direction					
		}	

		/*
		if(parameters.vDegree != undefined){
			//only raise the degree when initial geometry created.
			this.degreeRaise_V(); //attempt to raise the degree of the V direction to the specified value.			
		}	
		*/
		
		/////
		//W direction
		/////
		//If thickness parameters have been specified, make a solid.  If not, remain as a surface.
		//make through-thickness curves by offsetting and/or shrinking the surface.
		if(parameters.wDegree != undefined)this.wDegree = parameters.wDegree;
		if(parameters.centroidAtPercentChord !== undefined && parameters.centroidAtPercentChord>=0.)
			this.centroidAtPercentChord = parameters.centroidAtPercentChord;
			
		if(this.wDegree != undefined && this.wDegree>0){
			this.calculate_W_geometry(this.centroidAtPercentChord);
			//raise the degree of the model in the W direction to the specified value.		
			this.degreeRaise_W(); //only raise the degree when initial geometry created.
		}
		else if(this.wDegree == 0){
			//create a zero-thickness 'solid' to represent a surface.
			this.wDegree = 1;
			this.r = 1;
			this.W = [0,0,1,1];//linear knot vector;

			if(this.curves.length>1){
				//remove all previous geometry except for the original surface.
				this.curves.splice(1,this.curves.length); //.splice is a javascript array management method.
			}

			//make the 'inner' surface an exact copy of the surface curves.
			var wCurves = [];
			for(var i=0;i<this.curves[0].length;i++)wCurves.push(this.curves[0][i].clone());
			this.curves.push(wCurves);
		}
	}
	
	//copy constructor	
	clone(a){
		if(a === undefined)a = this;

		//var s = new NURBSSolid();
		var s = new this.constructor();

		s.curves = new Array(a.curves.length);
		for(var i=0;i<a.curves.length;i++){
			s.curves[i] = new Array(a.curves[0].length);
			for(var j=0;j<a.curves[i].length;j++){
				s.curves[i][j] = SurfaceCurve.prototype.clone(a.curves[i][j]);
			}
		}		

		s.q = a.q;
		s.r = a.r;	
		
		s.V = new Array(a.V.length);
		for(var i=0;i<s.V.length;i++)s.V[i]=a.V[i];
		
		s.W = new Array(a.W.length);
		for(var i=0;i<s.W.length;i++)s.W[i]=a.W[i];
		
		s.scale = a.scale;
		s.previousScale = a.previousScale;

		s.alignedAtPercentChord = a.alignedAtPercentChord;
		s.centroidAtPercentChord = a.centroidAtPercentChord;

		return s;
	}

	//scale the entire solid
	scaleSolid(scaleFactor){
		for(var i=0;i<this.curves.length;i++){
			for(var j=0;j<this.curves[i].length;j++){	
				this.curves[i][j].scaleCurve(scaleFactor/this.previousScale);
			}
		}	
		this.previousScale = scaleFactor;	
	}
	
	//rotate about a given point.
	rotateAboutPoint(point, rotationRadians){
		//rotate model about a 3D point. 
		var j;
		for(var i=0;i<this.curves.length;i++){
			for(j=0;j<this.curves[i].length;j++){
				this.curves[i][j].rotateAboutPoint(point, rotationRadians);
			}
		}	
	}
	
	//create geometry in the V direction.	
	calculate_knots(direction){
		//a routine to calculate the knots of the blade surface/thickness in the V or W direction
		//uses eq10.8 and eq9.8 from Piegl and Tiller.  This is the process known as "skinning".
		//input: direction = 0 for V; 1 for W.

		var vBar; //float[], our v_bar parameters
		var d; //float[], total chord lengths for each control point series along span
		
		var surface;
		if(direction==0)surface= this.curves[0]; //the array of surface curves in the V direction.
		else{
			//pack the surface array with the first row of curves in the W direction.
			surface = [];
			for(var i=0;i<this.curves.length;i++)surface.push(this.curves[i][0]);
		}
		
		/*
		//this.q = 3;	//degree of surface in V direction.  3 for cubic.
		var degree = 3; //degree of surface in V or W direction.  3 for cubic.
		if(surface[0].c.p<degree)degree=surface[0].c.p;
		if(surface.length==2)degree = 1; //straight ruled surface between two curves.
		*/
		//Default to the degree of the original curve.
		var degree = surface[0].c.p; //degree of surface in V or W direction.  3 for cubic.
		if(direction==0 && this.q>degree)degree=this.q;
		if(direction==1 && this.r>degree)degree=this.r;	
		
		//We need degree+1 curves for a surface of that degree.  Reduce degree if fewer curves.
		if(surface.length<degree+1)degree=surface.length-1;
		
		var m = surface.length - 1;   //high-index of control points in V or W direction
		var n = surface[0].c.n; //high index of control points in U direction	
		var k = m + degree + 1; //high index of knots in V or W direction

		//first find the total chord lengths for each series of control points along span
		var i, j, a, b;	
		d = new Array(n+1);
		for(i=0; i<=n; i++){
			d[i] = 0.;
			for(j=1; j<=m; j++){
				a = surface[j].c.P[i].p;
				b = surface[j-1].c.P[i].p;			
				d[i] += b.minus(a).NORM();
			}
		}

		//use eq10.8 to find the parameters v_k
		vBar = new Array(m+1);
		vBar[0] = 0.;
		vBar[m] = 1.;

		var temp;
		for(j=1; j<m; j++){
			temp = 0.;
			for(i=0; i<=n; i++){
				a = surface[j].c.P[i].p;
				b = surface[j-1].c.P[i].p;			
				
				temp += b.minus(a).NORM()/d[i];
			}

			//and calculate the parameter
			vBar[j] = vBar[j-1] + (1./(n+1))*temp;
		}
				
		//use eq9.8 to find the knots
		//this.V = new Array(k+1);
		var knotVector = new Array(k+1);
		for(j=0; j<=degree; j++)knotVector[j] = 0.;
		for(j=k-degree; j<=k; j++)knotVector[j] = 1.;

		for(j=1; j<=m-degree; j++){
			temp = 0.;
			for(i=j; i<=j+degree-1; i++)temp += vBar[i];
			
			knotVector[j+degree] = temp/degree;
		}
		
		if(direction == 0){
			//V direction
			this.q = degree;
			this.V = knotVector;
		}
		else{
			//W direction
			this.r = degree;
			this.W = knotVector;
		}
		
		return;	
	}
	
	//create geometry in the W direction.	
	calculate_W_geometry(centroidAtPercentChord){
		//input: an array of surface curves.
		//output:  a second array of surface curves aligned along the 'core' of the initial surface.
		//output:  also creates the W knot vector, and sets the degree of the inner W geometry, r, to 1 - linear.

		if(this.curves.length>1){
			//this is an update, so remove all previous geometry except for the original surface.
			this.curves.splice(1,this.curves.length); //.splice is a javascript array management method.
		}

		//determine whether we need to calculate the centroid
		var calculateCentroid = false; //it's been specified.
		if(centroidAtPercentChord === undefined || centroidAtPercentChord<0.)calculateCentroid = true;	
		
		//create internal geometry.
		var w_curves = new Array(this.curves[0].length);
		var j, centroid, numSteps, uStep1, uStep2, thickness, p1, p2, dist;
		for(var i=0;i<this.curves[0].length;i++){
			w_curves[i] = this.curves[0][i].clone();

			centroid = undefined;
			if(calculateCentroid){
				//we're not using a specified centroid.
				//Find the thickest point along the chord, and set the centroid in the middle of it.
				numSteps = 50;
				uStep1 = w_curves[i].u_LE/numSteps;
				uStep2 = (1.-w_curves[i].u_LE)/numSteps;
				thickness = 0.;
				for(j=1;j<numSteps;j++){
					p1 = w_curves[i].c.pointOnNURBSCurve(j*uStep1);
					p2 = w_curves[i].c.pointOnNURBSCurve(1.-(j*uStep2));
					dist = p1.minus(p2).NORM();
					if(dist>=thickness){
						//document.write(thickness+" | "+step+": "+i+","+j+": "+(uStep1*j)+","+(1.-(uStep2*j))+"<br>");										
						thickness=dist;
						centroid = p1.plus(p2.minus(p1).multiply(0.5));
						centroid.z = this.curves[0][i].c.P[0].p.z;
					}
				}
				
				//if the thickness in the middle of the curve is similar to the max thickness
				//place the curve there instead
				p1 = w_curves[i].c.pointOnNURBSCurve(w_curves[i].u_LE/2);
				p2 = w_curves[i].c.pointOnNURBSCurve(1.-(w_curves[i].u_LE/2));
				dist = p1.minus(p2).NORM();

				if(dist/thickness>0.9 && dist/thickness<1.1){
					thickness = dist;
					centroid = w_curves[i].c.P[0].p.plus(w_curves[i].LE.minus(w_curves[i].c.P[0].p).multiply(0.5));
					centroid.z = this.curves[0][i].c.P[0].p.z;
				}
			}

			//////////////////////
			//now we have a centroid, offset control points in reference to it.
			
			//if a point along the chord has been specified as the centroid, create the centroid from it.
			if(!centroid){
				centroid = this.curves[0][i].c.P[0].p.plus(this.curves[0][i].LE.minus(this.curves[0][i].c.P[0].p).multiply(centroidAtPercentChord));
				centroid.z = this.curves[0][i].c.P[0].p.z;
			}	

			//set all control points of the 'core' curve to this centroid.
			for(j=0;j<w_curves[i].c.P.length;j++){
				w_curves[i].c.P[j].p = centroid;
				w_curves[i].c.P[j].p.z = this.curves[0][i].c.P[0].p.z;
			}	
			w_curves[i].LE = centroid;
			w_curves[i].LE.z = this.curves[0][i].c.P[0].p.z;
		}
		this.curves.push(w_curves);
		
		//finally, create the W knot vector for this linear internal surface.
		this.calculate_knots(1); //1 == create surface/solid in W direction	
	}

	//raise the degree of the W direction.	
	degreeRaise_W(){
		if(this.wDegree<=this.r)return;
		//input:  a solid model with linear (degree 1) internal geometry, defined by a spanwise array of
		//surface curves, and a second spanwise array of 'core' curves degnerated to a point at the centroid 
		//of the surface.
		//output: a solid model with rebuilt curve arrays, featuring inner thickness curves to match the
		//specified degree of the W geometry, r.
		
		//make dummy curves for each control point in the W direction and elevate the degree of each.
		var i, j, k, m, new_r, new_W, newCurves;
		var numNewCurves = undefined, curve, x;
		for(i=0; i<this.curves[0].length;i++){ //for each surface curve along the span.
			for(j=0; j<this.curves[0][i].c.P.length;j++){ //for each control point in the curve.	
				//build a new W dummy curve, represented as a two-point line between a control point
				//on the surface and its corresponding point at the core.
				curve = new NURBSCurve();
				curve.n = 1;
				curve.p = 1;
				curve.P = new Array(2);
				curve.P[0] = this.curves[0][i].c.P[j]; //surface.
				curve.P[1] = this.curves[1][i].c.P[j]; //core.			
				curve.U = this.W;
				
				curve = curve.degreeElevateCurve(this.wDegree-1);
				for(k=0;k<curve.P.length;k++)curve.P[k].p.z = this.curves[0][i].c.P[0].p.z;

				//find how many new curves were created
				if(numNewCurves == undefined)numNewCurves = curve.P.length-2;
				//ensure that every new curve matches the initial degree elevation.
				else if(curve.P.length-2 != numNewCurves){
					//abort degree elevation.
					alert("ABORTED W DEGREE ELEVATION at curve[0]["+i+"].  Reverting to linear W geometry.");
					return;
				}

				//rebuild the curves array.
				if(i==0 && j==0){
					//the knot vector of the first elevated curve becomes the global W knot vector.
					new_W = [];
					for(k=0;k<curve.U.length;k++)new_W.push(curve.U[k]);
					//the degree of the first elevated curve becomes the global W degree.
					new_r = curve.p;
					//create a new curves array
					newCurves = [];
					//fill it with copies of the surface layer of curves, one for each control point in W curve.
					for(k=0;k<curve.P.length;k++){
						x = [];
						for(m=0;m<this.curves[0].length;m++)x.push(this.curves[0][m].clone());
						newCurves.push(x);
					}
				}
				//copy the positions of control points along the W direction.
				for(k=0;k<curve.P.length;k++){			
					newCurves[k][i].c.P[j] = curve.P[k].clone();
				}	
			}
		}
		
		//finalise the new geometry.
		this.W = new_W;
		this.r = new_r;
		//expand the curves array.  Don't remove the initial surface - it's hooked up to the GUI.
		for(i=1; i<newCurves.length; i++)this.curves.splice(i,1,newCurves[i]);
	}

	//raise the degree of the V direction.	
	degreeRaise_V(){
		if(this.q >= this.vDegree)return; //only attempt to elevate a lower degree geometry.
		//input:  a solid model with a variable number of spanwise curves.
		//output: a solid model with rebuilt curve arrays, featuring surface and inner curves to match the
		//specified degree of the V geometry, q.
		
		//make dummy curves for each control point in the V direction and elevate the degree of each.
		var i, j, k, m, new_q, new_V, newCurves;
		var numNewCurves = undefined, curve, x;
		for(i=0; i<this.curves.length;i++){ //for each curve along the W direction.
			for(j=0; j<this.curves[i][0].c.P.length;j++){ //for each control point in the first spanwise curve.
				//build a new V dummy curve, represented as a line of control points along the V direction.
				curve = new NURBSCurve();
				curve.n = this.curves[i].length-1;
				curve.p = this.q;
				curve.P = new Array(curve.n+1);
				for(k=0;k<this.curves[i].length;k++){
					curve.P[k] = this.curves[i][k].c.P[j].clone();
				}	
				curve.U = this.V;
				
				curve = curve.degreeElevateCurve(this.vDegree-this.q);

				//find how many new curves were created
				if(numNewCurves == undefined)numNewCurves = curve.P.length-this.curves[i].length;
				//ensure that every new curve matches the initial degree elevation.
				else if(curve.P.length-this.curves[i].length != numNewCurves){
					//abort degree elevation.
					alert("ABORTED V DEGREE ELEVATION at curve["+i+"]["+j+"].  Reverting to original V geometry.");
					return;
				}

				//rebuild the curves array.
				if(i==0 && j==0){
					//the knot vector of the first elevated curve becomes the global V knot vector.
					new_V = [];
					for(k=0;k<curve.U.length;k++)new_V.push(curve.U[k]);
					//the degree of the first elevated curve becomes the global V degree.
					new_q = curve.p;
					//create a new curves array
					newCurves = [];
					//fill it with copies of the first row of spanwise curves.
					for(m=0;m<this.curves.length;m++){
						x = [];
						for(k=0;k<curve.P.length;k++){
							x.push(this.curves[m][0].clone());
						}
						newCurves.push(x);	
					}
				}
				//copy the positions of control points along the V direction.
				for(k=0;k<curve.P.length;k++){			
					newCurves[i][k].c.P[j] = curve.P[k].clone();
				}	
			}
		}
		
		//we've only created new NURBSCurve objects along the span.  These need to be converted
		//to SurfaceCurve objects to make them fully compatible with everything else.
		//1. First and last original curves remain unaltered, since curve interpolates them.
		for(i=0;i<newCurves.length;i++){
			newCurves[i][0] = this.curves[i][0];
			newCurves[i][newCurves[i].length-1] = this.curves[i][this.curves[i].length-1];		
		}	
		
		//2. Update all SurfaceCurve parameters of curves away from end points
		var vec, baseChord;
		for(i=0; i<newCurves.length; i++){
			for(j=1; j<newCurves[0].length-1; j++){
				newCurves[i][j].vSpan = newCurves[i][j].c.P[0].p.z;
				newCurves[i][j].LE = newCurves[i][j].c.pointOnNURBSCurve(newCurves[i][j].u_LE);
				newCurves[i][j].LE = new Point_3D(newCurves[i][j].LE.x, newCurves[i][j].LE.y, newCurves[i][j].c.P[0].p.z);

				vec = newCurves[i][j].LE.minus(newCurves[i][j].c.P[0].p);			
				/*
				//instead of setting parameters like twist and scale to the actual values for this new curve in relation
				//to the base curve, redefine the curve as its own base curve.  This is the only thing that makes
				//sense for degree-elevated geometry - no base curve exists for new curves.
				newCurves[i][j].alphaDegrees = Math.atan(vec.y/vec.x)*180/Math.PI;
				*/
				
				baseChord = 1.;//define original as 1.0;
				
				newCurves[i][j].scale = vec.NORM()/baseChord;
			}
		}

		//finalise the new geometry.
		this.V = new_V;
		this.q = new_q;
		this.curves = newCurves;
	}

	//Compute 3D point on Bspline/NURBS surface.  NURBS just uses weighted points.	
	surfacePoint(u, v, index){
		//A4.3 in Piegl and Tiller.
		//input:  u, v (surface parameters in each direction.  float.)
		//output: S -> the point in 3D space representing S(u,v)
		if(index==undefined)index = 0;
		var numSlices = this.curves.length-1;
		if(numSlices<=0)numSlices=1;
		
		return this.solidPoint(u,v,index/numSlices); //find the point on the surface of the solid
	}

	//Compute 3D point on Bspline/NURBS surface.  NURBS just uses weighted points.	
	solidPoint(u, v, w){
		//A4.3 in Piegl and Tiller.
		//input:  u, v, w (surface parameters in each direction.  float.)
		//output: S -> the point in 3D space representing S(u,v,w)

		var i, j;
		
		//find the span and basis functions for the base curve at this u parameter
		////
		//u
		var uSpan = this.curves[0][0].c.findSpan(u);
		var Nu = this.curves[0][0].c.basisFuns(uSpan, u); //float[]
		////
		//v
		var m = this.curves[0].length-1; //high index of control points in the V direction
		var k = m+this.q+1;//high-index of knot vector in V direction

		//let's make a dummy bspline curve for calculations in the V direction
		var c = new NURBSCurve()
		c.initialise(m,this.q);
		for(i=0; i<=k; i++)c.U[i] = this.V[i];//make a copy of V knots

		var vSpan = c.findSpan(v); //float
		var Nv = c.basisFuns(vSpan, v); //float[]
		////
		//w
		var m_w = this.curves.length-1; //high index of control points in the W direction
		var k_w = m_w+this.r+1;//high-index of knot vector in W direction

		//a dummy bspline curve for calculations in the W direction
		var c_w = new NURBSCurve()
		c_w.initialise(m_w,this.r);
		for(i=0; i<=k_w; i++)c_w.U[i] = this.W[i];//make a copy of W knots

		var wSpan = c_w.findSpan(w); //float
		var Nw = c_w.basisFuns(wSpan, w); //float[]
		///////
		var p = this.curves[0][0].c.p;

		var temp = c_w.create2DArray(this.q+1, this.r+1);
		var iw, uInd, vInd, wInd;
		for(iw=0; iw<=this.r; iw++){
			wInd = wSpan - this.r + iw;	
			for(i=0; i<=this.q; i++){
				temp[i][iw] = new WeightedPoint_3D(0.,0.,0.,0.);

				vInd = vSpan - this.q + i;
				for(j=0; j<=p; j++){
					uInd = uSpan - p + j;		
					temp[i][iw] = temp[i][iw].plus(this.curves[wInd][vInd].c.P[uInd].pw().multiply(Nu[j]));//weighted point for NURBS.
				}
			}
		}	
		var s = new WeightedPoint_3D(0.,0.,0.,0.); //a point on our surface in (x,y,z,weight) space.
		for(iw=0; iw<=this.r; iw++){	
			for(i=0; i<=this.q; i++){
				s = s.plus(temp[i][iw].multiply(Nv[i]).multiply(Nw[iw]));
			}
		}	
		
		return s.divide(s.w).p;	//project the point down to x,y,z space and return it.	
	}

	//extract a 2D curve at a given v parameter
	extractSurfaceIsocurve(v){
		//input: v - parameter along V direction of surface.
		//output: NURBSCurve representing the slice of surface at the given v.
		//Eq 4.17 p135 in P&T.
		
		var i,j;
		
		var m = this.curves[0].length-1; //high index of control points in the V direction
		var k = m+this.q+1;//high-index of knot vector in V direction

		//make a dummy bspline curve for calculations in the V direction
		var c = new NURBSCurve()
		c.initialise(m,this.q);
		for(i=0; i<=k; i++)c.U[i] = this.V[i];//make a copy of V knots

		var vSpan = c.findSpan(v);
		var Nv = c.basisFuns(vSpan, v);

		//the isocurve will be similar to any of the original defining curves.
		var isoCurve = this.curves[0][0].c.clone();
		
		//initialise control points.
		for(i=0;i<isoCurve.P.length;i++)isoCurve.P[i] = new WeightedPoint_3D(0.,0.,0.,0.);
		
		var vInd;
		for(j=0; j<=this.q; j++){
			vInd = vSpan - this.q + j;
			for(i=0; i<isoCurve.P.length; i++){
				isoCurve.P[i] = isoCurve.P[i].plus(this.curves[0][vInd].c.P[i].pw().multiply(Nv[j]));
			}
		}

		//project down to (x,y,z) space.
		for(i=0;i<isoCurve.P.length;i++)isoCurve.P[i] = isoCurve.P[i].divide(isoCurve.P[i].w);
		
		return isoCurve;
	}

	//compute non-rational surface derivatives at the parameter values u and v	
	bsplineSurfaceDerivs(u, v, d){
		//A3.6 in Piegl and Tiller.  Pg 111.  'SurfaceDerivsAlg1'.
		//input:  u, v, d (highest order derivative to calculate)
		//output: SKL[][] (WeightedPoint_3D[][]), surface derivatives
		
		return this.bsplineSolidDerivs(u,v,0,d);
	}

	//compute non-rational surface derivatives at the parameter values u and v	
	bsplineSolidDerivs(u, v, w, d){
		//A3.6 in Piegl and Tiller.  Pg 111.  'SurfaceDerivsAlg1'.
		//input:  u, v, w, d (highest order derivative to calculate)
		//output: SKLA[][] (WeightedPoint_3D[][]), surface derivatives
		//K,L,A are the derivative orders of U,V,W directions, respectively.
		//So SKLA[0][0][0] refers to the zeroth derivative, which is the actual point at S(u,v,w)
		//SKLA[0][0][1] is dS/dw; SKLA[1][1][0] is d^2S/dudv, etc.
		
		var i, j, k;
				
		//ensure the number of derivatives specified is not more than the degree of each defining curve.
		var p = this.curves[0][0].c.p;
		var du=Math.min(d, p);
		var dv=Math.min(d, this.q);
		var dw=Math.min(d, this.r);	
		
		//make a dummy curve to characterise an isoline in the V direction
		var cv = new NURBSCurve()
		var m = this.curves[0].length-1; //high-index of control points in V direction
		cv.initialise(m,this.q);
		for(i=0; i<=this.V.length; i++)cv.U[i] = this.V[i];//make a copy of V knots

		//make a dummy curve to characterise an isoline in the W direction
		var cw = new NURBSCurve()
		var mw = this.curves.length-1; //high-index of control points in W direction
		cw.initialise(mw,this.r);
		for(i=0; i<=this.W.length; i++)cw.U[i] = this.W[i];//make a copy of W knots

		
		//find derivatives of the basis functions in each direction
		var uspan = this.curves[0][0].c.findSpan(u);
		var Nu = this.curves[0][0].c.dersBasisFuns(uspan, u, du); //float[][]
		var vspan = cv.findSpan(v);
		var Nv = cv.dersBasisFuns(vspan, v, dv);
		var wspan = cw.findSpan(w);
		var Nw = cw.dersBasisFuns(wspan, w, dw);
		
		//initialise our output derivatives array,  WeightedPoint_3D[][]
		var SKLA = this.create3DArray(du+1, dv+1, dw+1);
		
		for(i=0;i<=du;i++){
			for(j=0;j<=dv;j++){
				for(k=0;k<=dw;k++){
					SKLA[i][j][k] = new WeightedPoint_3D(0.,0.,0.,0.);
				}
			}
		}	
		
		var temp, su, sv, sw, dd, l, a, s;
		for(k=0; k<=du; k++){
			//var temp = new Array(this.q-1);
			temp = this.curves[0][0].c.create2DArray(this.q+1, this.r+1);

			for(sw=0; sw<=this.r; sw++){
				for(sv=0; sv<=this.q; sv++){
					temp[sv][sw] = new WeightedPoint_3D(0.,0.,0.,0.);
					
					for(su=0; su<=p;su++){
						temp[sv][sw] = temp[sv][sw].plus(this.curves[wspan-this.r+sw][vspan-this.q+sv].c.P[uspan-p+su].pw().multiply(Nu[k][su]));
					}
				}
			}	
				
			dd = Math.min(d-k, dv, dw);
			for(l=0;l<=dd; l++){
				for(a=0;a<=dd; a++){		
					for(s=0; s<=this.q; s++){
						for(sw=0; sw<=this.r; sw++){				
							SKLA[k][l][a] = SKLA[k][l][a].plus(temp[s][sw].multiply(Nv[l][s]).multiply(Nw[a][sw]));
						}	
					}
				}	
			}
		}
		
		//the elements of SKLA now contain A^(k,l,a) from equation 4.19 (as the 3D point), and W(u,v,w) as the weight.
		//SKLA[0][0][0] is just the 'zeroth derivative' in both directions, which is just the surface point S(u,v,w).
		//To get these derivatives in (x,y,z) space, divide by weight:  S(k,l,a).p/S(k,l,a).w
		//Remember:  this only gives derivatives of the non-rational Bspline surface.  
		//To get the surface derivaties, A4.4 needs to be used on SKLA.
		return SKLA;
	}

	//a utility function to create a 3D array.
	create3DArray(a, b, c){
		var i, j, k;
		var x = new Array(a);
		for(i=0;i<a;i++){
			x[i] = new Array(b);
			for(j=0;j<b;j++){
				x[i][j] = new Array(c);
				for(k=0;k<b;k++){
					x[i][j][k] = 0.;
				}
			}	
		}	
		return x
	}

	//find surface derivatives.
	rationalSurfaceDerivs(u, v, d){
		return this.rationalSolidDerivs(u,v,0,d);
	}

	//derivatives at any point in a NURBS solid
	rationalSolidDerivs(u, v, w, d){
		//A4.4 in Piegl and Tiller.
		//Input:  u,v,w,d -> bsplineSurfaceDerivs -> SKLA -> Aders and wders
		//output:  SKLA with rational derivatives
		var Aders = this.bsplineSolidDerivs(u, v, w, d);
		
		var p = this.curves[0][0].c.p;
		var du=Math.min(d, p);
		var dv=Math.min(d, this.q);	
		var dw=Math.min(d, this.r);

		var SKLA = this.create3DArray(du+1, dv+1, dw+1);

		var i, j, k;
		for(i=0;i<=du;i++){
			for(j=0;j<=dv;j++){
				for(k=0;k<=dw;k++){		
					SKLA[i][j][k] = new WeightedPoint_3D(0.,0.,0.,0.);
				}
			}
		}	

		var l, a, b, v1, v2, v3;
		for(k=0; k<=du; k++){
			for(l=0; l<=dv; l++){
				for(a=0; a<=dw; a++){		
					//Assemble equation 4.20 (p136) for S(u,v,w)			
					v1 = Aders[k][l][a];
					
					//u
					for(i=1; i<=k; i++){
						v1 = v1.minus(SKLA[k-i][l][a].multiply(Aders[i][0][0].w).multiply(this.binomial(k,i)));
					}
					//v
					for(i=1; i<=l; i++){
						v1 = v1.minus(SKLA[k][l-i][a].multiply(Aders[0][i][0].w).multiply(this.binomial(l,i)));
					}
					//w
					for(i=1; i<=a; i++){
						v1 = v1.minus(SKLA[k][l][a-i].multiply(Aders[0][0][i].w).multiply(this.binomial(a,i)));
					}
					//uv
					for(i=1; i<=k; i++){
						v2 = new WeightedPoint_3D(0.,0.,0.,0.);
						for(j=1; j<=l; j++){				
							v2 = v2.plus(SKLA[k-i][l-j][a].multiply(Aders[i][j][0].w).multiply(this.binomial(l,j)));
						}
						v1 = v1.minus(v2.multiply(this.binomial(k,i)));
					}
					//uw
					for(i=1; i<=k; i++){
						v2 = new WeightedPoint_3D(0.,0.,0.,0.);
						for(j=1; j<=a; j++){				
							v2 = v2.plus(SKLA[k-i][l][a-j].multiply(Aders[i][0][j].w).multiply(this.binomial(a,j)));
						}
						v1 = v1.minus(v2.multiply(this.binomial(k,i)));
					}
					//vw
					for(i=1; i<=l; i++){
						v2 = new WeightedPoint_3D(0.,0.,0.,0.);
						for(j=1; j<=a; j++){				
							v2 = v2.plus(SKLA[k][l-i][a-j].multiply(Aders[0][i][j].w).multiply(this.binomial(a,j)));
						}
						v1 = v1.minus(v2.multiply(this.binomial(l,i)));
					}
					
					//uvw
					for(i=1; i<=k; i++){
						v2 = new WeightedPoint_3D(0.,0.,0.,0.);
						for(j=1; j<=l; j++){
							v3 = new WeightedPoint_3D(0.,0.,0.,0.);
							for(b=1; b<=a; b++){				
								v3 = v3.plus(SKLA[k-i][l-j][a-b].multiply(Aders[i][j][b].w).multiply(this.binomial(a,b)));
							}
							v2 = v2.plus(v3.multiply(this.binomial(l,j)));
						}
						v1 = v1.minus(v2.multiply(this.binomial(k,i)));
					}
					
					//divide by weight to get derivative (k,l,a)
					SKLA[k][l][a] = v1.divide(Aders[0][0][0].w);
				}	
			}
		}
		return SKLA;
	}

	//Binomial lookup table via Pascal's triangle.
	binomial(n, k){
		//returns: int
		/*
		Instead of calculating the binomial, look it up from a table
		Here's Pascal's triangle for n<=4 from the wikipedia article 
		http://en.wikipedia.org/wiki/Binomial_coefficient:
		0:					1 								
		1:				1 		1 							
		2:			1 		2 		1 						
		3:		1 		3 		3 		1 					
		4:	1 		4 		6 		4 		1 				
		*/
		
		var bin = this.curves[0][0].c.create2DArray(5,5);
		
		bin[0][0]=1;
		
		bin[1][0]=1;
		bin[1][1]=1;
		
		bin[2][0]=1;
		bin[2][1]=2;
		bin[2][2]=1;
		
		bin[3][0]=1;
		bin[3][1]=3;
		bin[3][2]=3;
		bin[3][3]=1;
		
		bin[4][0]=1;
		bin[4][1]=4;
		bin[4][2]=6;
		bin[4][3]=4;
		bin[4][4]=1;
		
		return bin[n][k];
	}

	//area of NURBS outer surface via Gaussian quadrature.
	getSurfaceArea(u_start, u_end, v_start, v_end, u_steps_per_span, v_steps_per_span){
		//if no parameters have been specified, calculate area of entire surface
		if(!v_steps_per_span){
			u_start = this.curves[0][0].c.U[0];
			u_end = this.curves[0][0].c.U[this.curves[0][0].c.U.length-1];
			v_start = this.V[0];
			v_end = this.V[this.V.length-1];

			u_steps_per_span = 4;
			v_steps_per_span = 4;
		}

		var i, j;
		
		//assemble the vertices of our surface elements
		var uSpan = [], end, last_u, stepSize;
		uSpan.push(u_start);
		for(i = this.curves[0][0].c.p; i<=this.curves[0][0].c.n+1; i++){
			end = this.curves[0][0].c.U[i];
			if(u_end<end)end = u_end;
			if(this.curves[0][0].c.U[i]>u_start && end<=u_end){
				last_u = uSpan[uSpan.length-1];
				stepSize = (end - last_u)/u_steps_per_span;
				if(stepSize>0.){
					for(j=1;j<u_steps_per_span;j++)uSpan.push(last_u+(j*stepSize));
					uSpan.push(end);
				}	
			}	
		}

		var vSpan = [], last_v;
		vSpan.push(v_start);
		for(i = this.q; i<=this.curves[0].length; i++){	
			end = this.V[i];
			if(v_end<end)end = v_end;
			if(this.V[i]>v_start && end<=v_end){
				last_v = vSpan[vSpan.length-1];
				stepSize = (end - last_v)/v_steps_per_span;
				if(stepSize>0.){
					for(j=1;j<v_steps_per_span;j++)vSpan.push(last_v+(j*stepSize));
					vSpan.push(end);
				}	
			}	
		}
		
		//We're integrating {f(u) = |dS/du|.|dS/dv|} to find the area of the surface.
		//Perform Gauss quadrature to determine the area of differential elements, then sum.

		//there are n+p+2 knots in a curve, with 2p repeated end knots.
		//So the first knot in the span will be U[p], and the last U[n+p+1-p] = U[n+1]
		//Knots in the V direction are similar.
		//Each knot span represents one polynomial 'element', attached at the end by knots.
		//perform gauss quadrature to find the length of each of these elements in each direction
		
		//use the two-point gauss rule on each element.
		//employ a change of interval to shift integration from [-1,1] to [u_i, u_i+1]
		//https://en.wikipedia.org/wiki/Gaussian_quadrature
		var gauss = Math.sqrt(1./3.); // two-point gauss has two u, ±sqrt(1/3) on [-1,1]
		
		var surface = this.curves[0]; //the array of surface curves

		var sum = 0.;	
		var u_sum, u_a, u_b, v_a, v_b, u0, u1, v0, v1;
		var deriv00, deriv01, deriv10, deriv11, dA00, dA01, dA10, dA11;
		for(i=0; i<vSpan.length-1;i++){
			v_a = vSpan[i];
			v_b = vSpan[i+1];

			u_sum = 0.;		
			if(v_b>v_a){
				//two-point gauss rule.
				v0 = -gauss*(v_b-v_a)/2. + (v_a+v_b)/2.;
				v1 = gauss*(v_b-v_a)/2. + (v_a+v_b)/2.;	
				
				//assume a common knot vector for all curves in surface
				for(j=0; j<uSpan.length-1;j++){
					u_a = uSpan[j];
					u_b = uSpan[j+1];

					if(u_b>u_a){
						//two-point gauss rule.
						u0 = -gauss*(u_b-u_a)/2. + (u_a+u_b)/2.;
						u1 = gauss*(u_b-u_a)/2. + (u_a+u_b)/2.;	

						deriv00 = this.rationalSurfaceDerivs(u0,v0,1);
						deriv01 = this.rationalSurfaceDerivs(u0,v1,1);			
						deriv10 = this.rationalSurfaceDerivs(u1,v0,1);
						deriv11 = this.rationalSurfaceDerivs(u1,v1,1);
						
						dA00 = deriv00[0][1][0].NORM()*deriv00[1][0][0].NORM();
						dA01 = deriv01[0][1][0].NORM()*deriv01[1][0][0].NORM();			
						dA10 = deriv10[0][1][0].NORM()*deriv10[1][0][0].NORM();
						dA11 = deriv11[0][1][0].NORM()*deriv11[1][0][0].NORM();			
						
						u_sum += ((u_b-u_a)/2)*(dA00+dA01+dA10+dA11);
					}	
				}
				sum += ((v_b-v_a)/2)*u_sum;
			}	
		}
		
		return sum;
	}

	//volume of NURBS solid via Gaussian quadrature.
	getVolume(u_start, u_end, v_start, v_end, w_start, w_end, u_steps_per_span, v_steps_per_span, w_steps_per_span){
		if(!this.W || this.W.length<4)return 0.;
		
		//if no parameters have been specified, calculate volume of entire solid
		if(!w_steps_per_span){
			u_start = this.curves[0][0].c.U[0];
			u_end = this.curves[0][0].c.U[this.curves[0][0].c.U.length-1];
			v_start = this.V[0];
			v_end = this.V[this.V.length-1];
			w_start = this.W[0];
			w_end = this.W[this.W.length-1];
			
			u_steps_per_span = 4;
			v_steps_per_span = 4;
			w_steps_per_span = 4;		
		}
		
		//assemble the vertices of our volume elements
		//
		//there are n+p+2 knots in a curve, with 2p repeated end knots.
		//So the first knot in the span will be U[p], and the last U[n+p+1-p] = U[n+1]
		//Knots in the V and W directions are similar.
		//Each knot span represents one polynomial 'element', attached at the end by knots.
		//perform gauss quadrature to find the length of each of these elements in each direction
		
		var i, j, k, uSpan = [], end, last_u, stepSize;
		uSpan.push(u_start);
		for(i = this.curves[0][0].c.p; i<=this.curves[0][0].c.n+1; i++){
			end = this.curves[0][0].c.U[i];
			if(u_end<end)end = u_end;
			if(this.curves[0][0].c.U[i]>u_start && end<=u_end){
				last_u = uSpan[uSpan.length-1];
				stepSize = (end - last_u)/u_steps_per_span;
				if(stepSize>0.){
					for(j=1;j<u_steps_per_span;j++)uSpan.push(last_u+(j*stepSize));
					uSpan.push(end);
				}	
			}	
		}

		var vSpan = [], last_v;
		vSpan.push(v_start);
		for(i = this.q; i<=this.curves[0].length; i++){	
			end = this.V[i];
			if(v_end<end)end = v_end;
			if(this.V[i]>v_start && end<=v_end){
				last_v = vSpan[vSpan.length-1];
				stepSize = (end - last_v)/v_steps_per_span;
				if(stepSize>0.){
					for(j=1;j<v_steps_per_span;j++)vSpan.push(last_v+(j*stepSize));
					vSpan.push(end);
				}	
			}	
		}
		
		var wSpan = [], last_w;
		wSpan.push(w_start);
		for(i = this.r; i<=this.curves.length; i++){	
			end = this.W[i];
			if(w_end<end)end = w_end;
			if(this.W[i]>w_start && end<=w_end){
				last_w = wSpan[wSpan.length-1];
				stepSize = (end - last_w)/w_steps_per_span;
				if(stepSize>0.){
					for(j=1;j<w_steps_per_span;j++)wSpan.push(last_w+(j*stepSize));
					wSpan.push(end);
				}	
			}	
		}
		
		//We're integrating {f(u) = |dS/du|.|dS/dv|.|dS/dw|} to find the volume of the solid.
		//Perform Gauss quadrature to determine the volume of differential elements, then sum.
		var sum = 0.;
		
		//use the two-point gauss rule on each element.
		//employ a change of interval to shift integration from [-1,1] to [u_i, u_i+1]
		//https://en.wikipedia.org/wiki/Gaussian_quadrature
		var gauss = Math.sqrt(1./3.); // two-point gauss has two u, ±sqrt(1/3) on [-1,1]

		var u_a, u_b, v_a, v_b, w_a, w_b, u_sum, v_sum, u0, u1, v0, v1, w0, w1;
		var deriv000, deriv010, deriv100, deriv110, deriv001, deriv011, deriv101, deriv111;
		var dV000, dV010, dV100, dV110, dV001, dV011, dV101, dV111;	
		for(k=0; k<wSpan.length-1;k++){
			w_a = wSpan[k];
			w_b = wSpan[k+1];
			v_sum = 0.;
			if(w_b>w_a){
				//two-point gauss rule.
				w0 = -gauss*(w_b-w_a)/2. + (w_a+w_b)/2.;
				w1 = gauss*(w_b-w_a)/2. + (w_a+w_b)/2.;	
			
				for(i=0; i<vSpan.length-1;i++){
					v_a = vSpan[i];
					v_b = vSpan[i+1];

					u_sum = 0.;
					if(v_b>v_a){
						//two-point gauss rule.
						v0 = -gauss*(v_b-v_a)/2. + (v_a+v_b)/2.;
						v1 = gauss*(v_b-v_a)/2. + (v_a+v_b)/2.;	
						
						//assume a common knot vector for all curves in surface
						for(j=0; j<uSpan.length-1;j++){
							u_a = uSpan[j];
							u_b = uSpan[j+1];

							if(u_b>u_a){
								//two-point gauss rule.
								u0 = -gauss*(u_b-u_a)/2. + (u_a+u_b)/2.;
								u1 = gauss*(u_b-u_a)/2. + (u_a+u_b)/2.;	

								deriv000 = this.rationalSolidDerivs(u0,v0,w0,1);
								deriv010 = this.rationalSolidDerivs(u0,v1,w0,1);							
								deriv100 = this.rationalSolidDerivs(u1,v0,w0,1);
								deriv110 = this.rationalSolidDerivs(u1,v1,w0,1);
								
								deriv001 = this.rationalSolidDerivs(u0,v0,w1,1);
								deriv011 = this.rationalSolidDerivs(u0,v1,w1,1);							
								deriv101 = this.rationalSolidDerivs(u1,v0,w1,1);
								deriv111 = this.rationalSolidDerivs(u1,v1,w1,1);
								
								dV000 = deriv000[0][1][0].NORM()*deriv000[1][0][0].NORM()*deriv000[0][0][1].NORM();
								dV010 = deriv010[0][1][0].NORM()*deriv010[1][0][0].NORM()*deriv010[0][0][1].NORM();
								dV100 = deriv100[0][1][0].NORM()*deriv100[1][0][0].NORM()*deriv100[0][0][1].NORM();
								dV110 = deriv110[0][1][0].NORM()*deriv110[1][0][0].NORM()*deriv110[0][0][1].NORM();
								
								dV001 = deriv001[0][1][0].NORM()*deriv001[1][0][0].NORM()*deriv001[0][0][1].NORM();
								dV011 = deriv011[0][1][0].NORM()*deriv011[1][0][0].NORM()*deriv011[0][0][1].NORM();
								dV101 = deriv101[0][1][0].NORM()*deriv101[1][0][0].NORM()*deriv101[0][0][1].NORM();
								dV111 = deriv111[0][1][0].NORM()*deriv111[1][0][0].NORM()*deriv111[0][0][1].NORM();
								
								u_sum += ((u_b-u_a)/2)*(dV000+dV010+dV100+dV110+dV001+dV011+dV101+dV111);
							}	
						}
						v_sum += ((v_b-v_a)/2)*u_sum;
					}	
				}
				sum += ((w_b-w_a)/2)*v_sum;		
			}
		}	
		
		return sum;
	}

	
	output(){
		var output = [];
		for(var k=0; k<this.curves.length; k++){	
			output.push("<br>----------<br>");
			output.push("<br>----------<br><br>");		
			output.push("W thickness section:  "+k+"<br>");
		
			for(var j=0; j<this.curves[k].length; j++){
				output.push("<br>----------<br>");
				output.push("Curve:  "+j+"<br>");
				output.push("(n="+this.curves[k][j].c.n+", p="+this.curves[k][j].c.p+", scale="+this.curves[k][j].scale+")<br>");
				for(var i=0;i<=this.curves[k][j].c.n;i++)output.push("curves["+j+"].P["+i+"].x: "+this.curves[k][j].c.P[i].p.x+", .y: "+this.curves[k][j].c.P[i].p.y+", .z: "+this.curves[k][j].c.P[i].p.z+"<br>");
				for(var i=0;i<=this.curves[k][j].c.n+this.curves[k][j].c.p+1;i++)output.push("curves["+j+"].U["+i+"]:  "+this.curves[k][j].c.U[i]+"<br>");

				output.push("<br>");
			}

			output.push("<br>----------<br>");
			output.push("V direction: (m="+(this.curves[0].length-1)+", q="+this.q+")<br>");
			for(var i=0;i<this.V.length;i++)output.push("blade.V["+i+"]:  "+this.V[i]+"<br>");
			
			output.push("W direction: (m="+(this.curves.length-1)+", r="+this.r+")<br>");
			for(var i=0;i<this.W.length;i++)output.push("blade.W["+i+"]:  "+this.W[i]+"<br>");
			
			
			for(var i=0;i<output.length;i++)document.write(output[i]);
		}	
	}
}
//
//Web worker for area and volume calculations.
self.onmessage = function(e){
	if(e.data.solidModel){
		self.importScripts('Point_2D.js','Point_3D.js','WeightedPoint_2D.js','WeightedPoint_3D.js','NURBSCurve.js','SurfaceCurve.js');

		var solidModel = NURBSSolid.prototype.clone(e.data.solidModel);
		
		var area = solidModel.getSurfaceArea();
		
		var result = {"area": area, "volume": undefined};
		//console.log(result);
		
		self.postMessage({"result":result});		
		
		var volume = solidModel.getVolume();
		result = {"area": area, "volume": volume};
		
		//console.log(result);		
		
		self.postMessage({"result":result});
	}	
};

