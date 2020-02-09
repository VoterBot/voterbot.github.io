"use strict";
class NURBSCurve{
	constructor(){
		this.n = 0;//high-index of control point vector (n+1 control points)
		this.p = 0;//degree of curve

		this.P = [];//weighted control point vector WeightedPoint_2D[]
		this.U = [];//knot vector:  there are n+p+2 knots (m=n+p+1 is the high-index of knot vector) float[]
	}

	//methods
	//optional constructor.  Default object can be created just with new NURBSCurve();	
	initialise(n, p){
		this.n = n;
		this.p = p;

		this.P = new Array(this.n+1); //WeightedPoint_2D[]
		this.U = new Array(this.n+this.p+2); //float[]
	}

	//copy constructor
	clone(a){
		var c;
		if(a === undefined){
			c = new this.constructor(); //make the NURBSCurve class extensible.	
			a = this;
		}	
		else c = new NURBSCurve();
		c.n = a.n;
		c.p = a.p;

		c.P = new Array(a.n+1); //WeightedPoint_2D[]
		c.U = new Array(a.n+a.p+2); //float[]

		var i;	
		for(i=0;i<a.P.length;i++){
			if(a.P[i].p.z === undefined)c.P[i] = WeightedPoint_2D.prototype.clone(a.P[i]);
			else c.P[i] = WeightedPoint_3D.prototype.clone(a.P[i]);
		}
		for(i=0;i<a.U.length;i++){
			c.U[i] = a.U[i];
		}
		
		return c;
	}

	//load a bspline curve from file.  Set all control point weights to 1.0.
	createCurveFromFile(bsplineData){
		//split the data file into lines and store them in an array
		var data = bsplineData.split("\n");
		
		//first two lines contain n and p
		this.n = parseInt(data[0]);
		this.p = parseInt(data[1]);		

		//create new control point and knot arrays
		this.P = new Array(this.n+1);//control point vector <Point_2D>
		this.U = new Array(this.n+this.p+2);//knot vector <double> there are n+p+2 knots (m=n+p+1 is the high-index of knot vector)		

		var i;
		
		//load control points
		var coords;
		for(i=0;i<=this.n;i++){
			coords = data[i+2].split(','); //var coords
			//this.P[i] = new Point_2D(parseFloat(coords[0]),parseFloat(coords[1]));
			this.P[i] = new WeightedPoint_2D(parseFloat(coords[0]),parseFloat(coords[1]),1.);			
		}

		//load knots
		for(i=0;i<=this.n+this.p+1;i++)this.U[i] = parseFloat(data[i+this.n+3]);
	}  

	//Ex 4.2 from Piegl and Tiller - representation of a circle.
	setCircleCurve(){
		//first two lines contain n and p
		this.n = 8;
		this.p = 2;

		this.P = new Array(this.n+1); //WeightedPoint_2D[]
		this.U = new Array(this.n+this.p+2); //float[]
			
		//control points and weights for circle centred at 0,0, radius 1.
		this.P[0]=new WeightedPoint_2D(1.,0.,1.);
		this.P[1]=new WeightedPoint_2D(1.,1.,1.);
		this.P[2]=new WeightedPoint_2D(0.,1.,2.);
		this.P[3]=new WeightedPoint_2D(-1.,1.,1.);
		this.P[4]=new WeightedPoint_2D(-1.,0.,1.);
		this.P[5]=new WeightedPoint_2D(-1.,-1.,1.);
		this.P[6]=new WeightedPoint_2D(0.,-1.,2.);
		this.P[7]=new WeightedPoint_2D(1.,-1,1);
		this.P[8]=new WeightedPoint_2D(1.,0.,1.);
		
		//knots
		this.U[0] = 0.;
		this.U[1] = 0.;
		this.U[2] = 0.;
		this.U[3] = 0.25;
		this.U[4] = 0.25;
		this.U[5] = 0.5;	
		this.U[6] = 0.5;	
		this.U[7] = 0.75;
		this.U[8] = 0.75;		
		this.U[9] = 1.;
		this.U[10] = 1.;
		this.U[11] = 1.;
	}

	//determine the knot span index - which two knots in the U knot vector lie on either	
	findSpan(u){
		//ALGORITHM A2.1 in Piegel and Tiller - The NURBS Book
		//side of the given u parameter?
		//if(this.n+this.p+2 != this.U.length)return this.p;
		
		if(u == this.U[this.n+1])return this.n;//special case
		
		//do binary search
		var low=this.p;
		var high=this.n+1;
		var mid = Math.floor((low+high)/2);
		
		while(u< this.U[mid] || u>= this.U[mid+1]){
			if(u<this.U[mid])high=mid;
			else low=mid;
			
			mid = Math.floor((low+high)/2);
		}
		
		return mid;	
	}

	//Compute the nonvanishing basis functions - N
	basisFuns(span, u){
		//ALGORITHM A2.2 in Piegel and Tiller - The NURBS Book
		//returns float[]
		
		var N = new Array(this.p+1);
		var left = new Array(this.p+1);
		var right = new Array(this.p+1);
		
		N[0]=1.0;
		var r, saved, temp;
		for(var i=1;i<=this.p;i++){
			left[i] = u-this.U[span+1-i];
			right[i] = this.U[span+i]-u;
			saved = 0.0; //var saved
			for(r=0; r<i; r++){
				temp = N[r]/(right[r+1]+left[i-r]); //var temp
				N[r] = saved+right[r+1]*temp;
				saved = left[i-r]*temp;
			}
			N[i] = saved;	
		}
		
		return N;	
	}

	//create a 2D array
	create2DArray(a, b){
		var x = new Array(a), j;
		for(var i=0;i<a;i++){
			x[i] = new Array(b);
			for(j=0;j<b;j++)x[i][j]=0.;
		}	
		
		return x
	}

	//Compute non-zero basis functions and their derivatives.
	dersBasisFuns(span, u, upperDeriv){
		//ALGORITHM A2.2 in Piegel and Tiller - The NURBS Book
		//First section is A2.2 modified to store functions and knot differences.
		//returns float[][]

		//'span' is the index of the knot at the start of the section of the curve containing
		//the parameter 'u'.  'upperDeriv' is the maximum number of derivatives to calculate: d<=p
			
		var ndu = this.create2DArray(this.p+1, this.p+1);
		var ders = this.create2DArray(this.p+1, this.p+1);
		var a = this.create2DArray(2,this.p+1);
		var left = new Array(this.p+1);
		var right = new Array(this.p+1);

		var j,k;	
		
		ndu[0][0] = 1.0;
		var saved, r, temp;
		for(j=1; j<=this.p; j++){
			left[j] = u - this.U[span+1-j];
			right[j] = this.U[span+j] - u;
			saved = 0.0; //var saved

			for(r=0; r<j; r++){
				//lower triangle
				ndu[j][r] = right[r+1] + left[j-r];
				temp = ndu[r][j-1]/ndu[j][r]; //var temp
				//upper triangle
				ndu[r][j] = saved + right[r+1]*temp;
				saved = left[j-r]*temp;		
			}
			ndu[j][j] = saved;
		}

		//load the basis functions
		for(j=0; j<=this.p; j++)ders[0][j] = ndu[j][this.p];

		//this section computes the derivatives (eq2.9)
		var s1, s2, d, rk, pk, j1, j2;
		for(r=0; r<=this.p; r++){ //loop over function index
			s1 = 0; //var s1
			s2 = 1;	//var s2 - alternate rows in array a
			a[0][0] = 1.0;

			//loop to compute kth derivative
			for(k=1; k<=upperDeriv; k++){
				d = 0.0; //var d
				rk = r-k; //var rk
				pk = this.p-k; //var pk
				j1 = 0; //var j1
				j2 = 0; //var j2
				if(r>=k){
					a[s2][0] = a[s1][0]/ndu[pk+1][rk];
					d = a[s2][0]*ndu[rk][pk];
				}
				if(rk >= -1)j1 = 1;
				else j1 = -rk;
				if(r-1 <= pk)j2 = k-1;
				else j2 = this.p-r;

				for(j=j1; j<=j2; j++){
					a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pk+1][rk+j];
					d += a[s2][j]*ndu[rk+j][pk];
				}

				if(r <= pk){
					a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
					d += a[s2][k]*ndu[r][pk];
				}
			
				ders[k][r] = d;
				j=s1;
				s1=s2;
				s2=j; //switch rows
			}
		}		

		//multiply through by the correct factors (eq2.9)
		r = this.p;
		for(k=1; k<=upperDeriv; k++){
			for(j=0; j<=this.p; j++)ders[k][j]*=r;
			r *= (this.p-k);
		}

		//ders[k][j] is the kth derivative of the function N_i-p+j,p
		//where 0<=k<=upperDeriv and 0<=j<=p
				
		return ders;
	}

	//compute curve derivatives at the parameter value u
	rationalCurveDerivs(u, d){
		//A3.2 and A4.2
		//input:  n, p, C, u, d (highest order derivative to calculate)
		//output: Ck[] (Point_2D[])

		//ensure the number of derivatives specified is not more than the degree of the curve
		var du=this.p;
		if(d < this.p)du=d;
				
		//find derivatives of the basis functions
		var span = this.findSpan(u);
		var nders = this.dersBasisFuns(span, u, du); //float[][]
				
		//first, use the non-rational curve equations to find the derivatives of the numerator and the
		//denominator of equation 4.1 Piegl and Tiller
		var Aders = new Array(du+1); //Point_2D[]
		var wders = new Array(du+1); //float[]

		var i,j,k;
		
		var is2D = true;
		if(this.P[0].p.z !== undefined)is2D = false;
		
		//calculate derivatives of each non-rational component of the curve: A3.2
		var pw;
		for(k=0; k<=du; k++){
			if(is2D)Aders[k] = new Point_2D(0.,0.);
			else Aders[k] = new Point_3D(0.,0.,0.);		
			wders[k] = 0.;

			for(j=0; j<=this.p; j++){
				pw = this.P[span-this.p+j].pw(); //var pw - WeightedPoint_2D
				Aders[k] = Aders[k].plus(pw.p.multiply(nders[k][j]));
				wders[k] += pw.w*nders[k][j];
			}        
		}
				
		//now find rational derivatives: A4.2
		var Ck = new Array(du+1); //Point_2D[]
		var v;
		for(k=0; k<=du; k++){
			if(is2D)Ck[k] = new Point_2D(0.,0.);
			else Ck[k] = new Point_3D(0.,0.,0.);		
						
			v = Aders[k].clone(); //var v - Point_2D
						
			for(i=1; i<=k; i++){
				v = v.minus(Ck[k-i].multiply(this.binomial(k,i)*wders[i]));
			}
					  
			Ck[k] = v.divide(wders[0]);  
		}

		return Ck;
	}

	//compute a point on the NURBS curve: A4.1 Piegl and Tiller	
	pointOnNURBSCurve(u){
		var span = this.findSpan(u);
		var basisFunctions = this.basisFuns(span,u); //double[]

		var cw;
		if(this.P[0].p.z == undefined)cw = new WeightedPoint_2D(0.,0.,0.); //WeightedPoint_2D
		else cw = new WeightedPoint_3D(0.,0.,0.,0.); //WeightedPoint_3D	
		
		var pw;
		for(var i=0;i<=this.p;i++){
			//weight the control point
			pw = this.P[span-this.p+i].pw(); //var pw - WeightedPoint_2D
			cw = cw.plus(pw.multiply(basisFunctions[i]));
		}  
		
		var c = cw.divide(cw.w);//project 3D curve onto 2D space
		
		return c.p;
	}

	scale(xfac, yfac, zfac){
		for(var i=0;i<this.P.length;i++){
			this.P[i].p.x *= xfac;
			this.P[i].p.y *= yfac;
			if(zfac && this.P[i].p.z)this.P[i].p.z *= zfac;
		}  
	}

	translate(xMove, yMove, zMove){
		for(var i=0;i<this.P.length;i++){
			this.P[i].p.x += xMove;
			this.P[i].p.y += yMove;
			if(zMove && this.P[i].p.z)this.P[i].p.z += zMove;
		}  
	}

	//move the whole curve so that P[0] appears at the specified point.
	move(px, py, pz){
		if(pz && this.P[0].p.x == px && this.P[0].p.y == py && this.P[0].p.z == pz)return;	
		else if(pz===undefined && this.P[0].p.x == px && this.P[0].p.y == py)return;

		var xMove = px - this.P[0].p.x;
		var yMove = py - this.P[0].p.y;	
		var zMove = undefined;
		if(pz && this.P[0].p.z !== undefined)zMove =  pz - this.P[0].p.z;		
		
		for(var i=0;i<this.P.length;i++){
			this.P[i].p.x += xMove;
			this.P[i].p.y += yMove;
			if(zMove)this.P[i].p.z += zMove;		
		}  
	}

	//Project point2D q on to current curve to get Rk=C(uk)	
	projectPoint(curveOrPoint){
		// eq6.5 and eq6.4
		/////
		//Modified: input can now be a 2D point *or* a curve object.  Distinguish them by checking for .x variable.
		/////
		//q.output("q");
		//output ub and ek -> the parameter assoiated with this point and its error
		//returns: float[]
		//If input is a curve object, project its nearest point to the current curve.

		//partition the curve into num_segments segments.  The endpoints of each segment
		//are evaluated as potential "knot spans" for the parameter which is the
		//projection of Q[i].q to the curve.  Find the one that gives Rk=C(uk) closest to Q[i].q
		var i, j;
		
		var umin=this.U[0];
		var umax=this.U[this.U.length-1];
		
		var Rk = this.pointOnNURBSCurve(umin); //Initial approximation.
		
		//define storage for a candidate solution for point projection, or for an array of solutions for
		//curve intersection.
		var solution = [];
		
		var min_distance;
		var curveB = undefined;
		var q = undefined;
		//If the provided q is a point it will have q.x and q.y parameters.
		//If it doesn't have q.x, it's a curve.
		if(curveOrPoint.x === undefined){
			//curve.
			curveB = curveOrPoint;
		}
		else{
			//point.
			q = curveOrPoint;
			//min_distance = Rk.distance2D(q);
			min_distance = Rk.minus(q).NORM();

			//First solution.  Constitutes the only solution for point projection.		
			solution.push({
				u_A: umin,  
				error_A: undefined, 
				point_A: this.pointOnNURBSCurve(umin), //Point_2D
				//
				u_B:undefined, 
				error_B:undefined, 
				point_B:q.clone(),
			});
		}	

		//find a first approximation to the projected point / curve intersections.
		var tol=1.0e-5;	        
		var numSegments = 20; //segments per knot span.
		var distance;
		var a,b, segmentLength, Rk_minus, v1, v2, uq;
		
		//divide each unique knot span into linear segments.  Check segment end points.
		var spans = this.getKnotSpans();
		for(i=1;i<spans.length;i++){
			a = spans[i-1];
			b = spans[i];
			segmentLength = (b-a)/numSegments;
			Rk_minus = this.pointOnNURBSCurve(a);
			
			for(j=1;j<=numSegments;j++){
				ubar = a + j*segmentLength;
		
				if(ubar < umin)ubar = umin;
				if(ubar > umax)ubar = umax;
				Rk = this.pointOnNURBSCurve(ubar);
				
				if(curveB){
					//For curve intersections, find places where linear segments of curve A cross over curve B.		
					uq = curveB.projectPoint(Rk);

					v1 = Rk.minus(uq.point_A);
					v2 = Rk_minus.minus(uq.point_A);
					if(v1.plus(v2).NORM()<=v2.minus(v1).NORM()){
						var tmp = ubar - (segmentLength*v1.NORM()/(v1.NORM()+v2.NORM()));
						//the two ends of the linear segment of curveA straddle curveB.
						solution.push({
							u_A: tmp, 
							//u_A: ubar,
							error_A: undefined, 
							point_A: undefined,
							//
							u_B: undefined, 
							error_B: undefined, 
							point_B: undefined,
						});
					}	

					Rk_minus = Rk.clone();
				}	
				else{
					//For points, find the distance between the candidate point Rk and the field point q.
					//distance = Rk.distance2D(q);
					distance = Rk.minus(q).NORM();
				
					if(distance < min_distance){
						min_distance = distance;
						umin = ubar;
						solution[0].u_A = umin;
						solution[0].error_A = distance;
						solution[0].point_A = Rk.clone();				
					}
				}	
			}
		}	
		
		//if(curveB)console.log("intersection candidates: "+solution.length);
		
		//now that U_zero has been found, use newton iteration to find u_bar
		//within tolerance.
		var Ck; //Ck[] contains the zeroth, first and second derivatives of C(uk)		
		var numerator, denom1, denom2, ubar, iterations, ek;
		var maxIterations = 20;
		if(curveB)maxIterations = 4000;//curve intersection needs to be precise.
		for(i=0;i<solution.length;i++){
			umin = solution[i].u_A;
			if(umin < this.U[0])umin = this.U[0];
			if(umin > this.U[this.U.length-1])umin = this.U[this.U.length-1];

			for(iterations=0;iterations<maxIterations;iterations++){
				ubar = umin;
				//find the zeroth, first and second derives for eq6.3
				Ck = this.rationalCurveDerivs(ubar, 2);
				
				if(curveB){
					uq = curveB.projectPoint(Ck[0]); //zeroth derivative is the point on curveA.		
					q = uq.point_A;
				}	
				
				numerator = Ck[1].dot(Ck[0].minus(q));
				denom1 = 0.; //for linear curves.
				if(Ck[2])denom1 = Ck[2].dot(Ck[0].minus(q));
				denom2 = Math.pow(Ck[1].NORM(),2);
				
				umin = ubar - (numerator/(denom1+denom2));
				
				//now check for convergence using eq 6.4 criteria and break if met
				if(Ck[0].minus(q).NORM()<=tol)break; //criteria 1
				if(Math.abs(Ck[1].dot(Ck[0].minus(q)))/
					(Ck[1].NORM()*(Ck[0].minus(q)).NORM())<=tol)break; //criteria 2
				
				if(umin < this.U[0])umin = this.U[0];
				if(umin > this.U[this.U.length-1])umin = this.U[this.U.length-1];

				if(Ck[1].multiply(umin - ubar).NORM() <= tol)break;//criteria 3
			}//end of iteration loop.
			
			//iteration has convered to a u_bar parameter within tolerance, so update Q[i].ub and Q[i].ek

			//a small sanity check - sometimes floating point error results in small overruns at end points.
			if(Math.abs(umin - this.U[0]) < tol){
				umin = this.U[0];
				Ck[0] = this.P[0].p;
			}	
			if(Math.abs(umin - this.U[this.U.length-1]) < tol){
				umin = this.U[this.U.length-1];
				Ck[0] = this.P[this.P.length-1].p;
			}
			
			//ek = Ck[0].distance2D(q);
			ek = Ck[0].minus(q).NORM();
			
			//repack the solutions vector with the best solution, even if it is outside of tolerance.
			solution[i].u_A = umin;
			solution[i].error_A = ek;
			solution[i].point_A = Ck[0].clone();
			//find normal at this point.
			var normal = new Point_2D();
			normal.x = -Ck[1].y;
			normal.y = Ck[1].x;
			solution[i].normal_A = normal.divide(normal.NORM());
			
			//
			if(curveB){
				solution[i].u_B = uq.u_A;
				solution[i].error_B = uq.error_A;
				solution[i].point_B = uq.point_A.clone();
				solution[i].normal_B = uq.normal_A.clone();			
			}	
		}	

		if(curveB)return solution; //curve.
		return solution[0]; //point.
	}

	//modified version of projectPoint (above).  Takes a 2D point with an unknown x or y ordinate,
	//and finds the other ordinate.
	projectCoordinate(x, y){
		if(x !== undefined && y !== undefined)return point; //we need one unknown.

		var i, j;
		
		var umin=this.U[0];
		var umax=this.U[this.U.length-1];
		
		var Rk = this.pointOnNURBSCurve(umin); //Initial approximation.
		
		var min_distance;

		if(x !== undefined)min_distance = Math.abs(Rk.x - x);
		else min_distance = Math.abs(Rk.y - y);

		//First solution.
		var solution = {
			u_A: umin,  
			error_A: min_distance, 
			point_A: this.pointOnNURBSCurve(umin), //Point_2D
		};

		//find a first approximation to the projected point / curve intersections.
		var tol=1.0e-5;	        
		var numSegments = 20; //segments per knot span.
		var distance;
		var a, b, segmentLength, Rk_minus;
		
		//divide each unique knot span into linear segments.  Check segment end points.
		var spans = this.getKnotSpans();
		for(i=1;i<spans.length;i++){
			a = spans[i-1];
			b = spans[i];
			segmentLength = (b-a)/numSegments;
			Rk_minus = this.pointOnNURBSCurve(a);
			
			for(j=1;j<=numSegments;j++){
				ubar = a + j*segmentLength;
		
				if(ubar < umin)ubar = umin;
				if(ubar > umax)ubar = umax;
				Rk = this.pointOnNURBSCurve(ubar);
				
				//find the 1D distance between the candidate point Rk and the field point.
				if(x !== undefined)distance = Math.abs(Rk.x - x);
				else distance = Math.abs(Rk.y - y);
			
				if(distance < min_distance){
					min_distance = distance;
					umin = ubar;
					solution.u_A = umin;
					solution.error_A = distance;
					solution.point_A = Rk.clone();				
				}
			}
		}	
		
		//now that U_zero has been found, use newton iteration to find u_bar
		//within tolerance.
		var Ck; //Ck[] contains the zeroth, first and second derivatives of C(uk)		
		var numerator, denom1, denom2, ubar, iterations, ek;
		var maxIterations = 20;

		umin = solution.u_A;
		if(umin < this.U[0])umin = this.U[0];
		if(umin > this.U[this.U.length-1])umin = this.U[this.U.length-1];

		for(iterations=0;iterations<maxIterations;iterations++){
			ubar = umin;
			//find the zeroth, first and second derives for eq6.3
			Ck = this.rationalCurveDerivs(ubar, 2);

			if(x !== undefined){
				numerator = Ck[1].x*(Ck[0].x - x);
				denom1 = 0.; //for linear curves.
				if(Ck[2])denom1 = Ck[2].x*(Ck[0].x - x);
				denom2 = Ck[1].x*Ck[1].x;
			}	
			else{
				numerator = Ck[1].y*(Ck[0].y - y);
				denom1 = 0.; //for linear curves.
				if(Ck[2])denom1 = Ck[2].y*(Ck[0].y - y);
				denom2 = Ck[1].y*Ck[1].y;
			}	

			umin = ubar - (numerator/(denom1+denom2));
			
			//now check for convergence using eq 6.4 criteria and break if met
			if(x !== undefined){
				if(Math.abs(Ck[0].x - x)<=tol)break; //criteria 1
				if(Math.abs(Ck[1].x*(Ck[0].x - x))/Math.abs(Ck[1].x*(Ck[0].x - x))<=tol)break; //criteria 2
			}
			else{
				if(Math.abs(Ck[0].y - y)<=tol)break; //criteria 1
				if(Math.abs(Ck[1].y*(Ck[0].y - y))/Math.abs(Ck[1].y*(Ck[0].y - y))<=tol)break; //criteria 2
			}
			
			if(umin < this.U[0])umin = this.U[0];
			if(umin > this.U[this.U.length-1])umin = this.U[this.U.length-1];

			if(x !== undefined){
				if(Math.abs(Ck[1].x*(umin - ubar)) <= tol)break;//criteria 3
			}
			else{
				if(Math.abs(Ck[1].y*(umin - ubar)) <= tol)break;//criteria 3
			}
				
		}//end of iteration loop.

		//iteration has convered to a u_bar parameter within tolerance, so update Q[i].ub and Q[i].ek

		//a small sanity check - sometimes floating point error results in small overruns at end points.
		if(Math.abs(umin-this.U[0])<tol){
			umin = this.U[0];
			Ck[0] = this.P[0].p;
		}	
		if(Math.abs(umin-this.U[this.U.length-1])<tol){
			umin = this.U[this.U.length-1];
			Ck[0] = this.P[this.P.length-1].p;
		}

		if(x !== undefined)ek = Math.abs(Ck[0].x - x);
		else ek = Math.abs(Ck[0].y - y);	
		
		//repack the solutions vector with the best solution, even if it is outside of tolerance.
		solution.u_A = umin;
		solution.error_A = ek;
		solution.point_A = Ck[0].clone();

		return solution;
	}

	//Pascal's triangle for binomial lookup.
	binomial(n, k){
		//returns: int
		/*
		Instead of calculating the binomial, look it up from a table
		Here's Pascal's triangle for n<=8 from the wikipedia article 
		http://en.wikipedia.org/wiki/Binomial_coefficient:
		0: 									1 								
		1: 								1 		1 							
		2: 							1 		2 		1 						
		3: 						1 		3 		3 		1 					
		4: 					1 		4 		6 		4 		1 				
		5: 				1 		5 		10 		10 		5 		1 			
		6: 			1 		6 		15 		20 		15 		6 		1 		
		7: 		1  		7  		21 		35 		35 		21 		7  		1  	
		8: 	1  		8  		28 		56 		70 		56 		28 		8  		1 
		*/
		
		var bin = this.create2DArray(9,9);
		
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
		
		bin[5][0]=1;
		bin[5][1]=5;
		bin[5][2]=10;
		bin[5][3]=10;
		bin[5][4]=5;
		bin[5][5]=1;

		bin[6][0]=1;
		bin[6][1]=6;
		bin[6][2]=15;
		bin[6][3]=20;
		bin[6][4]=15;
		bin[6][5]=6;
		bin[6][6]=1;	

		bin[7][0]=1;
		bin[7][1]=7;
		bin[7][2]=21;
		bin[7][3]=35;
		bin[7][4]=35;
		bin[7][5]=21;
		bin[7][6]=7;
		bin[7][7]=1;
		
		bin[8][0]=1;
		bin[8][1]=8;
		bin[8][2]=28;
		bin[8][3]=56;
		bin[8][4]=70;
		bin[8][5]=56;
		bin[8][6]=28;
		bin[8][7]=8;
		bin[8][8]=1;	

		return bin[n][k];
	}
	
	//two routines for transforming a NURBSCurve from local to screen coordinates.
	getScreenTransform(w, h, fac){
		//find the max width and height
		var minX = 1.E6;
		var minY = 1.E6;
		var maxX = -1.E6;
		var maxY = -1.E6;

		var x,y;
		for(var i=0;i<this.P.length;i++){
			x = this.P[i].p.x;
			y = this.P[i].p.y;
			
			if(x>maxX)maxX=x;
			if(y>maxY)maxY=y;
			if(x<minX)minX=x;
			if(y<minY)minY=y;
		}
		
		var curveWidth = maxX-minX;
		var curveHeight = maxY-minY;

		//scale by the smallest dimension, so everything
		//fits in the screen whilst maintaining the aspect ratio
		//of the original curve
		var scaleX = fac*w/curveWidth;
		var scaleY = fac*h/curveHeight;
		var scaler = Math.min(scaleX, scaleY);
		
		//line up the curve with the centre of the screen
		var xOffset = (w-curveWidth*scaler)/2.;
		var yOffset = (h-curveHeight*scaler)/2.;

		return [xOffset, yOffset, minX, minY, maxX, maxY, scaler];
	}

	toScreenCoordinates(screenTransform){
		var xOffset = screenTransform[0];
		var yOffset = screenTransform[1];
		var minX = screenTransform[2];
		var minY = screenTransform[3];
		var maxX = screenTransform[4];
		var maxY = screenTransform[5];		
		var scaler = screenTransform[6];

		var x = this.clone();		

		//Screen_Y0 is at the top of the screen.		
		//flip about the xaxis through middle of the Y range
		var midY = minY + (maxY - minY)/2;		

		var tmp = maxY;
		maxY = midY-minY;
		minY = midY-tmp;
		
		var i;
		for(i=0;i<x.P.length;i++){
			x.P[i].p.y = midY - x.P[i].p.y;
		}  
		
		//scale and translate
		for(i=0;i<x.P.length;i++){
			x.P[i].p.x = (x.P[i].p.x-minX)*scaler + xOffset;
			x.P[i].p.y = (x.P[i].p.y-minY)*scaler + yOffset;
		}  
		
		return x;
	}

	getLength(){
		//We're integrating {f(u) = |dC/du|} to find the length of the curve.
		//Perform Gauss quadrature to determine length of curve segments, then sum.
		var sum = 0.;
		//there are n+p+2 knots in a curve, with 2p repeated end knots.
		//So the first knot in the span will be U[p], and the last U[n+p+1-p] = U[n+1]
		//each knot span represents one polynomial 'element', attached at the end by knots.
		//perform gauss quadrature to find the length of each of these elements.
		
		//use the two-point gauss rule on each element.
		//employ a change of interval to shift integration from [-1,1] to [u_i, u_i+1]
		//https://en.wikipedia.org/wiki/Gaussian_quadrature
		var gauss = Math.sqrt(1./3.); // two-point gauss has two u, ±sqrt(1/3) on [-1,1]

		var a, b, u1, u2, deriv1, deriv2;
		for(var i = this.p; i<this.n+1; i++){		
			a = this.U[i]; //var a
			b = this.U[i+1]; //var b
			//document.write(i+": "+a+", "+b+"<br>");
			//don't add the contributions of zero-length knot spans
			if(b>a){
				//two-point gauss rule.
				u1 = -gauss*(b-a)/2. + (a+b)/2.; //var u1
				u2 = gauss*(b-a)/2. + (a+b)/2.;	 //var u2
				//document.write(i+": "+u1+", "+u2+"<br>");
				deriv1 = this.rationalCurveDerivs(u1,2)[1].NORM(); //var deriv1
				deriv2 = this.rationalCurveDerivs(u2,2)[1].NORM(); //var deriv2
				sum += ((b-a)/2)*(deriv1+deriv2);
			}	
		}
		
		return sum;
	}

	//Elevate the degree of a NURBS curve t times.
	degreeElevateCurve(t){
		//A5.9 in Piegl and Tiller
		//Returns:  a new NURBSCurve of elevated degree.
		if(t<1)return this.clone();
		
		//var c = new NURBSCurve();
		//var c = new this.constructor(); //make the NURBSCurve class extensible.
		var c = this.clone();
		
		var m = this.n + this.p + 1;
		c.p = this.p+t;
		var ph2 = Math.floor(c.p/2.);
		
		//compute Bezier degree elevation coefficients
		var bezalfs = this.create2DArray(this.p+t+1, this.p+1);
		bezalfs[0][0] = 1.;
		bezalfs[c.p][this.p] = 1.;

		var i, j, k;
		
		var inv, mpi;
		for(i=1; i<=ph2; i++){
			inv = 1./this.binomial(c.p, i); //var inv
			mpi = Math.min(this.p, i); //var mpi
			
			for(j=Math.max(0,i-t); j<=mpi; j++){
				bezalfs[i][j] = inv*this.binomial(this.p,j)*this.binomial(t, i-j);
			}
		}
		
		for(i=ph2+1; i<=c.p-1; i++){
			mpi = Math.min(this.p,i); //var mpi
			for(j=Math.max(0,i-t); j<=mpi; j++){
				bezalfs[i][j] = bezalfs[c.p-i][this.p-j];
			}
		}

		var mh = c.p;
		var kind = c.p+1;
		
		var r = -1;
		var a = this.p;
		var b = this.p+1;
		var cind = 1;
		var ua = this.U[0];
		
		c.P = [];  //new control point array
		c.P.push(this.P[0].pw()); //all calculations take place in (x,y,w) space

		c.U = [];  //new knot vector	
		for(i=0; i<=c.p; i++){
			c.U.push(ua);
		}	
		
		//is this a 2D or 3D curve?
		var newWeightedPoint;
		if(this.P[0].p.z === undefined)newWeightedPoint = new WeightedPoint_2D(0.,0.,0.);
		else newWeightedPoint = new WeightedPoint_3D(0.,0.,0.,0.);
		
		//initialise first Bezier segment
		var bpts = new Array(this.p+1);
		for(i=0;i<=this.p;i++){
			bpts[i] = this.P[i].pw();
		}	

		var nextBpts = new Array(this.p - 1);
		for(i=0;i<nextBpts.length;i++){
			//nextBpts[i] = new WeightedPoint_2D(0.,0.,0.);
			nextBpts[i] = newWeightedPoint.clone();
		}
		
		var ebpts = new Array(this.p+t+1);			
		for(i=0;i<ebpts.length;i++){
			//ebpts[i] = new WeightedPoint_2D(0.,0.,0.);
			ebpts[i] = newWeightedPoint.clone();
		}
		
		var tmp, mul, ub, oldr, lbz, rbz, numer, alfs, save, s;
		var first, last, den, bet, tr, kj, alf, gam;
		while(b < m){
			tmp = b; //var tmp
			while(b < m && this.U[b] == this.U[b+1])b++;
			
			mul = b-tmp+1; //var mul
			mh += mul+t;
			ub = this.U[b]; //var ub
			oldr = r; //var oldr
			r = this.p-mul;
			
			//insert knot u(b) r times
			lbz = 1; //var lbz
			if(oldr > 0)lbz = Math.floor((oldr+2)/2.);
			rbz = c.p; //var rbz
			if(r > 0)rbz = Math.ceil(c.p-(r+1)/2.); //GODDAMMIT PIEGL AND TILLER.  Ceiling, not floor.
				
			if(r > 0){
				//insert knots to get Bezier segment
				numer = ub-ua; //var numer
				alfs = new Array(this.p); //var alfs
				for(k=this.p; k>mul; k--){
					alfs[k-mul-1] = numer/(this.U[a+k] - ua);
				}	
				
				for(j=1; j<=r; j++){
					save = r-j; //var save
					s = mul+j; //var s
					
					for(k=this.p; k>=s; k--){
						bpts[k] = bpts[k].multiply(alfs[k-s]).plus(bpts[k-1].multiply(1.-alfs[k-s]));
					}
					nextBpts[save] = bpts[this.p].clone();				
				}
			}//end of 'insert knot'	
			
			//degree elevate Bezier
			for(i=lbz; i<=c.p; i++){
				//Only points lbz, ... , c.p are used below.		
				//ebpts[i] = new WeightedPoint_2D(0.,0.,0.);		
				ebpts[i] = newWeightedPoint.clone();
				mpi = Math.min(this.p,i); //var mpi
				for(j=Math.max(0,i-t); j<=mpi; j++){
					ebpts[i] = ebpts[i].plus(bpts[j].multiply(bezalfs[i][j]));
				}
			}//end of degree elevating Bezier

			if(oldr > 1){
				//must remove knot u=U[a] oldr times.
				first = kind-2; //var first
				last = kind; //var last
				den = ub-ua; //var den
				bet = (ub-c.U[kind-1])/den; //var bet
				
				for(tr=1; tr<oldr; tr++){
					//knot removal loop
					i = first; //var i
					j = last; //var j
					kj = j-kind+1; //var kj
					
					while(j-i > tr){
						//loop and compute the new control points for one removal step.
						if(i < cind){
							alf = (ub-c.U[i])/(ua-c.U[i]); //var alf
							c.P[i] = c.P[i].multiply(alf).plus(c.P[i-1].multiply(1.-alf));
						}
						if(j >= lbz){
							if(j-tr <= kind-c.p+oldr){
								gam = (ub-c.U[j-tr])/den; //var gam
								ebpts[kj] = ebpts[kj].multiply(gam).plus(ebpts[kj+1].multiply(1.-gam));
							}
							else ebpts[kj] = ebpts[kj].multiply(bet).plus(ebpts[kj+1].multiply(1.-bet));
						}
						
						i++;
						j--;
						kj--;
					}
					first--;
					last++;
				}
			}//end of removing knot, u=U[a].
			
			//load knot ua
			if(a != this.p){
				for(i=0;i<c.p-oldr; i++){
					c.U.push(ua);
					kind++;
				}
			}
			
			//load control points into c.P
			for(j=lbz; j<=rbz; j++){
				c.P.push(ebpts[j].clone());
				cind++;
			}
			
			//set up for next pass through loop.
			if(b < m){
				for(j=0;j<r;j++)bpts[j] = nextBpts[j].clone();
				for(j=r;j<=this.p;j++)bpts[j] = this.P[b-this.p+j].pw();
				a=b;
				b++;
				ua = ub;
			}
			else{
				//end knot.
				for(i=0; i<=c.p; i++)c.U.push(ub);
			}
		}//end of while-loop (b < m).
		
		//finalise curve.
		c.n = mh-c.p-1;

		for(i=0;i<c.P.length;i++)c.P[i].p = c.P[i].p.divide(c.P[i].w); //project down to (x,y) space.
		
		//c.output("c")	
		
		return c;
	}

	//Insert the knot vector X into existing knot vector.
	refineKnotVector(X){
		//A5.4 Piegl and Tiller.
		//This is the process known as 'knot refinement'.
		//Input: X - knot vector to insert.
		//Output:  newly calculated curve.
		
		//var c = new NURBSCurve();
		//var c = new this.constructor(); //make the NURBSCurve class extensible.
		var c = this.clone();
		
		var r = X.length - 1;
		var m = this.n + this.p + 1;
		var a = this.findSpan(X[0]);
		var b = this.findSpan(X[r]);
		b++;

		var i, j, k, l;
		
		c.P = new Array(this.n + r + 2);
		for(j=0;j<=a-this.p;j++)c.P[j] = this.P[j].pw(); //all calculations take place in (x,y,w) space.
		for(j=b-1;j<=this.n;j++)c.P[j+r+1] = this.P[j].pw();
		
		c.U = new Array(this.U.length + X.length); //the merged knot vector.
		for(j=0;j<=a;j++)c.U[j] = this.U[j];
		for(j=b+this.p; j<=m; j++)c.U[j+r+1] = this.U[j];
		
		i = b+this.p-1; //var i
		k = b+this.p+r; //var k
		
		var ind, alfa;
		for(j=r; j>=0; j--){
			while(X[j] <= this.U[i] && i>a){
				c.P[k-this.p-1] = this.P[i-this.p-1].pw();
				c.U[k] = this.U[i];
				k--;
				i--;
			}
			
			c.P[k-this.p-1] = c.P[k-this.p].clone();
			for(l=1;l<=this.p; l++){
				ind = k-this.p+l; //var ind
				alfa = c.U[k+l] - X[j]; //var alfa
				if(Math.abs(alfa) == 0.)c.P[ind-1] = c.P[ind].clone();
				else{
					alfa = alfa/(c.U[k+l] - this.U[i-this.p+l]);
					c.P[ind-1] = c.P[ind-1].multiply(alfa).plus(c.P[ind].multiply(1.-alfa));
				}
			}
			
			c.U[k] = X[j];
			k--;
		}
		
		//finalise curve
		c.n = c.P.length - 1;
		c.p = this.p;
		
		//project down into (x,y) space.
		for(i=0;i<c.P.length;i++)c.P[i].p = c.P[i].p.divide(c.P[i].w);
		
		return c;
	}

	//insert p+1 knots to break the curve at each intersection location.	
	split(intersections){
		var i, j;
		var tol=1.E-3;	
		var newKnots=[], u, numInsertions, multiplicity;
		var breaks = [];
		for(i=0;i<intersections.length;i++){					
			u = intersections[i];
			numInsertions = this.p+1;
			//check for existing knots at this parameter.
			for(j=this.p;j<=this.n+1;j++){
				multiplicity = this.findKnotMultiplicity(j);
				if(Math.abs(u - this.U[j])<tol){
					u = this.U[j];
					//count how many others knots exist at this multiplicity.
					multiplicity = this.findKnotMultiplicity(j);
					numInsertions -= multiplicity[1];
					j = multiplicity[0]+1;
				}
			}
			if(u-this.U[0]>tol && this.U[this.U.length-1]-u>tol && numInsertions>=0){
				breaks.push(u);
				for(j=0;j<numInsertions;j++){
					newKnots.push(u);
				}
			}	
		}
		
		//if(newKnots.length==0)return undefined;
		if(breaks.length==0)return undefined;
		//console.log(this.name+"\n----\n"+this.toString()+"\n------------\n");
		//for(var i=0;i<newKnots.length;i++)console.log("newKnots["+i+"]="+newKnots[i])
		
		//create a new curve by inserting this new knot vector.
		var refinedCurve;
		if(newKnots.length>0)refinedCurve = this.refineKnotVector(newKnots);
		else refinedCurve = this;
		
		breaks.push(refinedCurve.U[refinedCurve.U.length-1]);//make the final knot a break.
		
		//console.log(this.name+"\n------------\n"+refinedCurve.toString());	

		//break the curve at p+1 repeated knots, creating an array of subcurves.
		var newCurves = [];
		var first=0, index, first_u, last_u;
		for(i=0;i<breaks.length;i++){
			newCurves.push(this.clone());
			//newCurves.push(new this.constructor());
			//newCurves[i].p = this.p
			
			newCurves[i].U = [];
			index = first;
			while(true){
				u = refinedCurve.U[index];
				if(u<=breaks[i])newCurves[i].U.push(u);
				else break;
				index++;				
			}
			index--;
			
			//normalise knot vector.
			first_u = newCurves[i].U[0];
			last_u = newCurves[i].U[newCurves[i].U.length-1] - first_u;
			for(j=0;j<newCurves[i].U.length;j++){
				newCurves[i].U[j] = (newCurves[i].U[j]-first_u)/last_u;
			}	
			
			newCurves[i].P = [];
			for(j=first; j<index-refinedCurve.p; j++){
				newCurves[i].P.push(refinedCurve.P[j].clone());
			}	
			newCurves[i].n = newCurves[i].P.length-1;
			
			//set up for next curve
			first = index-this.p;
		}	
		
		return newCurves;
	}

	removeAllRemovableKnots(tolerance){
		var newCurve = this.clone();
		
		if(!tolerance)tolerance = 1;
		
		var i = newCurve.p+1; //first removable knot.
		var multiplicity, num, x;
		while(i<=newCurve.n){
			multiplicity = newCurve.findKnotMultiplicity(i); //var multiplicity
			i = multiplicity[0];
			num = multiplicity[1]; //var num
			x = newCurve.removeKnot(i, num, tolerance); //var x
			
			i = i-x[1]+1;
			newCurve = x[0].clone();
			
			//if(i>=newCurve.U.length-newCurve.p-1)break;
		}

		if(newCurve.n != newCurve.P.length-1)
			console.log("ERROR - NURBSCurve.removeAllRemovableKnots - c.n+1="+(newCurve.n+1)+" c.P.length="+newCurve.P.length+"\n");
			
		if(newCurve.n+newCurve.p+2 != newCurve.U.length)
			console.log("ERROR - NURBSCurve.removeAllRemovableKnots - n+p+2="+(newCurve.n+newCurve.p+2)+" c.U.length="+newCurve.U.length+"\n");
			
		
		return newCurve;
	}

	//find s, the multiplicity of the knot.	
	findKnotMultiplicity(r){
		//first, make sure r is the last knot with this value.
		var i;
		var r1 = r;
		for(i=r1; i<this.U.length; i++){	
			if(this.U[i] == this.U[r])r=i;
			else if(this.U[i]>this.U[r])break;
		}
		
		//count repeated knots at this value
		var s = 1; //a multiplicity of one means knot is unique.
		for(i=r-1; i>=0; i--){
			if(this.U[i] == this.U[r])s++;
			else if(this.U[i]<this.U[r])break;
		}
		
		return [r,s];
	}

	//return the number of distinct knot spans in this knot vector.	
	getKnotSpans(){
		var spans=[], multiplicity, u;
		for(var i=this.p;i<=this.n+1;i++){
			u = this.findKnotMultiplicity(i);
			i=u[0];
			spans.push(this.U[i]);
		}
		return spans;
	}

	//Attempt to remove knot u == U[r] num times.	
	removeKnot(r,num, TOL){
		//A5.8 in Piegl and Tiller.
		//Input: r, num
		//Output:  new curve; and t = number of successful knot removals.
		
		var i, j, k;
		
		//If tolerance not specified in argument list, calculate it.
		if(!TOL){
			//find variables for tolerance calculation.  eq(5.30).
			var d = 1.E-4; //arbitrary bound on deviation.
			//find the minimum weight on the original curve.
			var wmin = 1E6;
			for(i=0;i<this.P.length;i++){
				if(this.P[i].w<wmin)wmin = this.P[i].w;
			}
			
			//find the maximum (2D) distance between control points.  Will need to be 3D when 3D points used.
			var maxDist = 0., dist;
			for(i=1;i<this.P.length;i++){
				dist = this.P[i].p.minus(this.P[0].p).NORM(); //var dist
				if(dist>maxDist)maxDist=dist;
			}
			//assemble the tolerance.
			TOL = (d*wmin)/(1.+maxDist);
		}	
		
		//find s, the multiplicity of the knot.
		var multiplicity = this.findKnotMultiplicity(r);
		var r = multiplicity[0];	
		var s = multiplicity[1];
		
		var m = this.n+this.p+1;
		var ord = this.p+1;
		var fout = Math.floor((2*r-s-this.p)/2); //first control point out.
		var last = r-s;
		var first = r-this.p;
		
		var c = this.clone(); //our modified curve.
		for(i=0;i<c.P.length;i++){
			c.P[i] = c.P[i].pw(); ////all calculations occur in (x,y,w) space.
		}	
		var u = this.U[r];
		//This loop is equation 5.28
		var t, off, temp, ii, jj, remflag, badWeight, alfi, alfj, dist, x;
		if(num>s)num=s;//can only remove knots that exist.
		for(t=0; t<num; t++){
			off = first - 1; //var off - diff in index between temp and P
			temp = new Array(2*this.p+1);// var temp
			temp[0] = c.P[off].clone(); 
			temp[last+1-off] = c.P[last+1].clone();

			i = first; //var i
			j = last; //var j
			ii = 1; //var ii
			jj = last-off; //var jj
			remflag = 0; //var remflag
			
			//compute new control points for one removal step
			badWeight = false; //var badWeight
			while(j-i > t){
				alfi = (u - c.U[i])/(c.U[i+ord+t]-c.U[i]); //var alfi
				alfj = (u - c.U[j-t])/(c.U[j+ord]-c.U[j-t]); //var alfj
				temp[ii] = c.P[i].minus(temp[ii-1].multiply(1.-alfi)).divide(alfi);
				temp[jj] = c.P[j].minus(temp[jj+1].multiply(alfj)).divide(1.-alfj);
				
				//Disallow weights that are less than or equal to zero.
				if(temp[ii].w<=0.)badWeight=true;
				else if(temp[jj].w<=0.)badWeight=true;
				
				i++;
				ii++;
				j--;
				jj--;
			}//end while-loop.
			
			//check if knot is removable.
			if(!badWeight && j-i < t){
				dist = temp[ii-1].distance3D(temp[jj+1]); //var dist
				if(dist <= TOL)remflag = 1;
			}
			else if(!badWeight){
				alfi = (u - c.U[i])/(c.U[i+ord+t] - c.U[i]); //var alfi
				x = temp[ii+t+1].multiply(alfi).plus(temp[ii-1].multiply(1.-alfi)); //var x
				
				dist = c.P[i].distance3D(x); //var dist
				if(dist <= TOL)remflag = 1;
			}
			
			if(remflag == 0)break; //cannot remove any more knots.
			else{
				//Successful removal.  Save new control points.
				i = first;
				j = last;
				
				while(j-i > t){
					c.P[i] = temp[i-off].clone();
					c.P[j] = temp[j-off].clone();
					i++;
					j--;
				}
			}
			first--;
			last++;
		}//end of for-loop.
		if(t == 0){
			return [this, 0];
		}	
		
		for(k=r+1; k<=m; k++)c.U[k-t] = c.U[k]; //shift knots.
		c.U.splice(c.U.length-t-1, t); //remove the last 't' knots	
		
		j=fout; //var j
		i=j; //var i - P[j] through P[i] will be overwritten.
		
		for(k=1; k<t; k++){
			if(k%2==1)i++;
			else j--;
		}
		
		var tmp = [];
		for(k=0; k<j; k++)tmp.push(c.P[k]);
		for(k=i+1; k<=c.n; k++){
			tmp.push(c.P[k].clone()); //shift
			j++;
		}
		c.P = tmp;

		//finalise curve
		c.n = c.P.length - 1;
		c.p = this.p;

		//project down into (x,y) space.
		for(i=0;i<c.P.length;i++){
			c.P[i].p = c.P[i].p.divide(c.P[i].w);
		}	
		
		return [c, t];
	}

	//reverse the knot and control point vectors.
	reverse(){
		var i;
		var x = this.clone();
		
		//control points
		for(i=0;i<this.P.length;i++)x.P[i] = this.P[this.P.length-1-i].clone();
		
		//knots
		var maxU = this.U[this.U.length-1];
		for(i=0;i<this.U.length;i++)x.U[i] = maxU - this.U[this.U.length-1-i];
		
		return x;
	}

	//Arbitrary quadratic circular arcs in NURBS.
	makeNURBSCircle(centre, radius, startAngle, endAngle){
		//A7.1 in Piegl and Tiller.
		var i, j, index;
		
		//set defaults - return a circle with radius==1.0 centred at origin.
		if(centre===undefined){
			centre = new Point_2D(0.,0.);
			radius = 1.;
			startAngle = 0.;
			endAngle = 2*Math.PI;
		}
		
		if(endAngle < startAngle)endAngle = 2.*Math.PI + endAngle;
		var theta = endAngle - startAngle;
		
		var narcs = 4; //number of arcs.  4 for a full circle.
		if(theta <=Math.PI/2.)narcs = 1;
		else if(theta<=Math.PI)narcs = 2;
		else if(theta<=(3.*Math.PI/2.))narcs = 3;
		
		var dtheta = theta/narcs;
		this.n = 2*narcs; //n+1 control points.
		this.p = 2;//quadratic.
		
		var w1 = Math.cos(dtheta/2.); //dtheta/2 is base angle.

		//initialise starting vectors.
		var P0 = new Point_2D(0.,0.), T0 = new Point_2D(0.,0.);		
		P0.x = radius*Math.cos(startAngle);
		P0.y = radius*Math.sin(startAngle);
		P0 = centre.plus(P0);
		
		T0.x = -Math.sin(startAngle);
		T0.y = Math.cos(startAngle);
		
		this.P = new Array(this.n+1);
		this.P[0] = new WeightedPoint_2D(P0.x, P0.y, 1.);
		
		var angle=startAngle, P2 = new Point_2D(0.,0.), T2 = new Point_2D(0.,0.);
		index = 0;
		for(i=1; i<=narcs; i++){
			angle += dtheta;
			P2.x = radius*Math.cos(angle);
			P2.y = radius*Math.sin(angle);
			P2 = centre.plus(P2);
			this.P[index+2] = new WeightedPoint_2D(P2.x, P2.y, 1.);

			T2.x = -Math.sin(angle);
			T2.y = Math.cos(angle)
			
			//find the intersection of these two tangent lines.
			this.P[index+1] = this.findVectorIntersection(P0,T0,P2,T2);
			this.P[index+1].w = w1;
			
			index+=2;
			if(i<narcs){
				P0 = P2.clone();
				T0 = T2.clone();
			}
		}
		
		j = 2*narcs+1; //load the knot vector.
		this.U = new Array(this.n+this.p+2);
		for(i=0; i<3; i++){
			this.U[i] = 0.;
			this.U[i+j] = 1.;
		}
		
		switch(narcs){
			case 1: break;
			case 2:
				this.U[3] = this.U[4] = 0.5;
				break;
			case 3:
				this.U[3] = this.U[4] = 1./3.;
				this.U[5] = this.U[6] = 2./3.;
				break;
			case 4:
				this.U[3] = this.U[4] = 0.25;
				this.U[5] = this.U[6] = 0.5;
				this.U[7] = this.U[8] = 0.75;
				break;
		}
		
		//fix rounding errors.
		for(var i=0;i<this.P.length;i++){
			this.P[i].p.x = this.P[i].p.x.round(8);
			this.P[i].p.y = this.P[i].p.y.round(8);		
		}
	}

	//find the intersection between two vectors using Cramer's rule.
	findVectorIntersection(P0,T0,P2,T2){
		//see: http://math.stackexchange.com/questions/406864/intersection-of-two-lines-in-vector-form
		
		//l1=(x1,y1)+a(u1,v1); l2=(x2,y2)+b(u2,v2)

		var a = T0.NORM();
		var v1 = T0.divide(a);
		var b = T2.NORM();
		var v2 = T2.divide(b);

		//intersection when l1===l2	
		//|v1.x -v2.x| · |a| = |P2.x-P0.x|
		//|v1.y -v2.y|   |b|   |P2.y-P0.y|
		
		//Solve with Cramer's rule: https://en.wikipedia.org/wiki/Cramer's_rule
		var c1 = P2.x - P0.x;
		var c2 = P2.y - P0.y;
		var a1 = v1.x;
		var a2 = v1.y;
		var b1 = -v2.x;
		var b2 = -v2.y;
		
		var a = (c1*b2 - b1*c2)/(a1*b2 - b1*a2);
		//var b = (a1*c2 - c1*a2)/(a1*b2 - b1*a2);	
		
		var intersection = P0.plus(v1.multiply(a));
		
		return new WeightedPoint_2D(intersection.x, intersection.y, 1);
	}
	
	//create a simple curve from control points and the simplest possible knot vector.
	createSimpleCurve(requestedDegree){
		var x = this.clone();
		if(x.n<1)return x;
		
		var i;
		
		if(x.p==0)x.p=1;
		if(x.p>requestedDegree)x.p=requestedDegree;
		
		x.U = new Array(x.n+x.p+2);
		for(i=0;i<=x.p;i++){
			x.U[i] = 0.;
			x.U[x.U.length-i-1] = 1.;	
		}
		if(x.n>=1){
			var step = 1./x.n;
			var index = 1;
			for(i=x.p+1;i<=x.n;i++){
				x.U[i]=step*index;
				index++;
			}	
		}
		
		if(requestedDegree<=1 || x.n<requestedDegree){
			x.p=1;
			return x;
		}	

		x = x.degreeElevateCurve(requestedDegree-x.p);
		
		return x;
	}	

	toString(){
		var i;
		var out = ("n = "+this.n+"\n");
		out += ("p = "+this.p+"\n");
		var z;
		for(i=0;i<this.P.length;i++){
			z = "";
			if(this.P[i].p.z !== undefined)z = (", "+this.P[i].p.z);
			out += ("P["+i+"] = ("+this.P[i].p.x+", "+this.P[i].p.y + z + ", "+this.P[i].w+")\n");
		}
		for(i=0;i<this.U.length;i++){
			out += ("U["+i+"] = "+this.U[i]+"\n");
		}
		
		return out;
	}	

	output(){
		document.write(this.toString().replace(/\n/g,"<br>"));
	}	
}
//
Number.prototype.round = function(places) {
  var result = +(Math.round(this + "e+" + places)  + "e-" + places);
  if(isNaN(result))return 0.;
  return result;
}