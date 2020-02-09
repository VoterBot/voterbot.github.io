"use strict";
class SurfaceCurve{
	constructor(){
		this.baseCurve = new NURBSCurve(); //the NURBSCurve object on which his curve is based. chord=1. alpha=0.
		this.c = undefined; //the working curve.  Scaled, rotated, translated.
		this.u_LE; //the parameter at the leading edge of the curve (ie: the point furthest from the TE).
		this.LE_;  //Point_2D representation of the leading edge
		this.scale; //a scaling factor
		this.vSpan; //the location of this curve along the span of the surface.
		
		this.showControlPoints = true;
		this.showCurve = true;
	}	

	//methods
	initialiseFromFile(dataFile, scale_, alphaDegrees_, vSpan_, degree_){
		this.baseCurve.createCurveFromFile(dataFile);
		this.initialise(scale_, alphaDegrees_, vSpan_, degree_);
	}
	
	initialiseFromCurve(baseCurve_, scale_, alphaDegrees_, vSpan_, degree_){
		this.baseCurve = baseCurve_.clone();
		this.initialise(scale_, alphaDegrees_, vSpan_, degree_);
	}
	
	initialise(scale_, alphaDegrees_, vSpan_, degree_){
		//make a copy of the base curve and use it as the working curve
		this.c = this.baseCurve.clone();

		//defaults.
		if(scale_ === undefined)scale_ = 1.;
		if(alphaDegrees_ === undefined)alphaDegrees_ = 0.;
		if(vSpan_ === undefined)vSpan_ === 0.;
		if(degree_ === undefined)degree_ === this.c.p;

		this.findLeadingEdge(); //find the leading edge (this.u_LE and this.LE) for the working curve.
		this.scale = scale_;
		this.vSpan = vSpan_;
		
		//raise the degree of the base curve if it has been requested.
		//this should only need to be done once for each original curve.
		if(degree_ && degree_>this.c.p){
			this.c = this.c.degreeElevateCurve(degree_-this.c.p);
		}	
		
		//translate the curve so that the leading edge is at (0,0)	
		//var a = this.moveOriginTo(this.LE);
		var a = this.moveOriginToPercentChord(0.);//0% chord is the leading edge.
		
		//define the distance between the TE and LE, the chord length, as 1.0
		var chordLength = this.c.P[0].p.minus(this.LE).NORM();

		//scale the curve with the given scaling factor	
		var scalingFactor = this.scale/chordLength;
		for(var i=0;i<this.c.P.length;i++)
			this.c.P[i].p = this.c.P[i].p.multiply(scalingFactor);
		this.LE = this.LE.multiply(scalingFactor);	

		//convert the 2D curve into a 2.5D working curve (ie 2D with a Z span)			
		for(var i=0;i<this.c.P.length;i++){
			this.c.P[i] = new WeightedPoint_3D(this.c.P[i].p.x, this.c.P[i].p.y, this.vSpan, this.c.P[i].w);
			//this.c.P[i].output(i);		
		}
		this.LE = new Point_3D(this.LE.x, this.LE.y, this.vSpan);

		this.rotateAboutPoint(this.LE, [alphaDegrees_*Math.PI/180.,0.,0.]); //radians.  [xy,xz,yz]rotate about Z axis in XY plane.
		var alpha = -(Math.atan(this.c.P[0].p.y/this.c.P[0].p.x)*180./Math.PI);	
	}
	
	//copy constructor
	clone(a){
		if(a === undefined)a = this;
		
		var x = new SurfaceCurve();
		x.baseCurve = NURBSCurve.prototype.clone(a.baseCurve);
		x.c = NURBSCurve.prototype.clone(a.c);	
		x.u_LE = a.u_LE;
		if(a.LE.z === undefined)x.LE = new Point_2D(a.LE.x, a.LE.y);
		else x.LE = new Point_3D(a.LE.x, a.LE.y, a.LE.z);
		x.scale = a.scale;
		x.vSpan = a.vSpan;
		x.direction = a.direction;

		return x;
	}
	
	//rotationRadians == [xy, xz, yz];
	rotateAboutPoint(point, rotationRadians){
		var tol = 1.E-6;
		var rotations = [undefined, undefined, undefined];	
		for(var i=0;i<3;i++){
			if(Math.abs(rotationRadians[i])>tol){
				rotations[i] = {
					costh: Math.cos(rotationRadians[i]),
					sinth: Math.sin(rotationRadians[i]),
				}	
			}	
		}	
		
		//https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations	
		var p1, alpha;
		for(var i=0;i<this.c.P.length;i++){	
			if(rotations[0]){
				p1 = this.c.P[i].p.minus(point);		
				//xy plane = rotate about z axis.				
				this.c.P[i].p.x = (p1.x*rotations[0].costh - p1.y*rotations[0].sinth) + point.x;
				this.c.P[i].p.y = (p1.x*rotations[0].sinth + p1.y*rotations[0].costh) + point.y;
			}	
			if(rotations[1]){
				p1 = this.c.P[i].p.minus(point);		
				//xz plane = rotate about y axis.								
				this.c.P[i].p.x = (p1.x*rotations[1].costh + p1.z*rotations[1].sinth) + point.x;
				this.c.P[i].p.z = (-p1.x*rotations[1].sinth + p1.z*rotations[1].costh) + point.z;
			}	
			if(rotations[2]){
				p1 = this.c.P[i].p.minus(point);		
				//yz plane = rotate about x axis.								
				this.c.P[i].p.y = (p1.y*rotations[2].costh - p1.z*rotations[2].sinth) + point.y;
				this.c.P[i].p.z = (p1.y*rotations[2].sinth + p1.z*rotations[2].costh) + point.z;
			}
		}
	}

	//find the point on the curve furthest away from the trailing edge.
	findLeadingEdge(curve){
		//Use Newton-Raphson to iterate: https://en.wikipedia.org/wiki/Newton's_method.
		//Equations 6.3 and 6.4 in Piegl and Tiller, as used in NURBSCurve.ProjectPoint()

		if(curve === undefined)curve = this.c;
		
		//Chord is a line between the trailing edge (C(0.)) and the furthest point on the curve away from it.
		//Point furthest away is Leading Edge (C(u_LE)), and is defined by parameter u_LE.
		//u_LE is found when f(u) = C'(u)·(C(u) - C(0.)) = 0.
		var ubar = 0.5; //initial guess.
		
		var Ck, numerator, denom1, denom2, maxIterations=200, umin=ubar, tol=1.E-6;
		for(var i=0;i<maxIterations;i++){
			ubar = umin;
			Ck = curve.rationalCurveDerivs(ubar, 2);
			
			numerator = Ck[1].dot(Ck[0].minus(curve.P[0].p));
			denom1 = 0.; //for linear curves.
			if(Ck[2])denom1 = Ck[2].dot(Ck[0].minus(curve.P[0].p));
			denom2 = Math.pow(Ck[1].NORM(),2);
			
			umin = ubar - (numerator/(denom1+denom2));
			
			//now check for convergence using eq 6.4 criteria and break if met
			if((Ck[0].minus(curve.P[0].p)).NORM()<=tol)break; //criteria 1
			if(Math.abs(Ck[1].dot(Ck[0].minus(curve.P[0].p)))/
				((Ck[1]).NORM()*(Ck[0].minus(curve.P[0].p)).NORM())<=tol)break; //criteria 2
			
			if(umin < curve.U[0])umin = curve.U[0];
			if(umin > curve.U[curve.U.length-1])umin = curve.U[curve.U.length-1];

			if((Ck[1]).multiply(umin - ubar).NORM() <= tol)break;//criteria 3
		}

		if(this.c){
			this.u_LE = umin;
			this.LE = Ck[0].clone();
		}
		else return {u_LE: umin, point:Ck[0].clone()};
	}

	//scale	
	scaleCurve(scalingFactor){
		this.scale *= scalingFactor;
		for(var i=0;i<this.c.P.length;i++){
			//this.c.P[i].p = this.c.P[i].p.minus(this.LE).multiply(scalingFactor).plus(this.LE);
			this.c.P[i].p = this.c.P[i].p.multiply(scalingFactor);
		}
		this.LE = this.LE.multiply(scalingFactor);
	}
	
	moveOriginToPercentChord(percentChord){
		var newOrigin = this.c.P[0].p.plus(this.LE.minus(this.c.P[0].p).multiply(1.-percentChord));
		var a = this.moveOriginTo(newOrigin);
	}

	//translate the curve so that the specified point is at (0,0)	
	moveOriginTo(point){
		var z;
		if(this.c.P[0].p.z)z = this.c.P[0].p.z;
		
		for(var i=0;i<this.c.P.length;i++){
			this.c.P[i].p = this.c.P[i].p.minus(point);
			if(z)this.c.P[i].p.z = z;
		}	
		this.LE = this.LE.minus(point);
		if(z)this.LE.z = z;	
	}
	
	flipX(){
		//The original curves are specified with a trailing edge at (1,0) and leading edge at (0,0), 
		//with control points numbered anti-clockwise from trailing edge.
		//Flip the curve so that TE is(0,0) and LE (1,0).  This also flips control point numbering
		//to a **clockwise** order, which screws up the definition of the outside of the curve - flipping
		//control points also turns the curve inside-out.
		
		var pointA = this.baseCurve.P[0].p.clone();
		var pointB;
		if(this.c)pointB = this.c.P[0].p.clone();	
		for(var i=0; i<this.baseCurve.P.length;i++){
			//set trailing edge as origin.
			this.baseCurve.P[i].p = this.baseCurve.P[i].p.minus(pointA);
			if(this.c)this.c.P[i].p = this.c.P[i].p.minus(pointB);
			
			//flip X coordinate.
			this.baseCurve.P[i].p.x *= -1;
			if(this.c)this.c.P[i].p.x *= -1;
		}
	}
}	