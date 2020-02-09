"use strict";
//extend the NURBSSolid object class
class BsplineBlade extends NURBSSolid{
	constructor(){
		//call the NURBSSolid constructor
		super();

		//aerofoils.
		this.baseFoil;//a template curve (chord = 1.0, leading edge at origin, twist=0.0) - SurfaceCurve.
		this.hubFoil; //SurfaceCurve

		//parameters
		this.numFoils; //number of Bspline aerofoils along span.
		this.rBlade;  //the length of the blade; - float
		this.pitchAngle; //float
		this.transitionPercent; //float
		this.modelTransition; //boolean
		//
		this.hasTransition; //boolean - does the currently-defined blade have a transition piece?
		this.currentTransitionPercent; //float
		//
		this.chordAndTwist;
		//
		this.mesh; //panel method mesh object.

		this.outwardNormalDirection = 1; //1==base curves defined anti-clockwise; -1==clockwise.
	}

	//methods
	clone(a){
		//make a clone of the NURBSSolid object's elements.
		
		if(a === undefined)a = this;
		else a.constructor = BsplineBlade;
		
		var x = super.clone(a);
		
		x.baseFoil = SurfaceCurve.prototype.clone(a.baseFoil);
		x.hubFoil = SurfaceCurve.prototype.clone(a.hubFoil);	
		//
		x.rBlade = a.rBlade;
		x.numFoils = a.numFoils;
		x.pitchAngle = a.pitchAngle;
		x.transitionPercent = a.transitionPercent;
		x.modelTransition = a.modelTransition;
			
		//
		if(a.mesh)x.mesh = BladeMesh.prototype.clone(a.mesh);
		x.outwardNormalDirection = a.outwardNormalDirection;

		//
		if(a.chordAndTwist){
			x.chordAndTwist = {};
			x.chordAndTwist.template = [];
			var point;
			for(var i=0; i<a.chordAndTwist.template.length; i++){
				point = WeightedPoint_3D.prototype.clone(a.chordAndTwist.template[i]);
				if(a.chordAndTwist.template[i].d)point.d = WeightedPoint_3D.prototype.clone(a.chordAndTwist.template[i].d);
				x.chordAndTwist.template.push(point);
			}	
			x.chordAndTwist.bspline = LeastSquaresCurve.prototype.clone(a.chordAndTwist.bspline);	
		}	
		
		return x;
	}
	
	createBlade(settings){
		this.setBladeSettings(settings);
		this.initialiseModel({isReset: true});
	}
	
	setBladeSettings(settings){
		this.numFoils = settings.numFoils;
		this.pitchAngle=settings.pitchAngle;
		this.rBlade=settings.rBlade;
		this.transitionPercent=settings.transitionPercent;
		this.modelTransition=settings.modelTransition;
		
		if(settings.baseFoil && settings.hubFoil){
			this.baseFoil = settings.baseFoil.clone();
			this.hubFoil = settings.hubFoil.clone();
		
			//ensure the aerofoils have a)origin at the trailing edge; and b)a positive X axis.
			this.baseFoil.flipX();
			this.hubFoil.flipX();
			this.outwardNormalDirection = -1; //mirroring about Y changes control point numbering to clockwise.
		}

		this.chordAndTwist = undefined;
		this.curves = undefined;
	}
	
	initialiseModel(chordAndTwistParameters){
		//console.log("BsplineBlade.js: initialiseModel()");
		
		//If chordAndTwistParameters have been passed as an argument, update blade's definition.
		if(chordAndTwistParameters)this.getChordAndTwistDistribution(chordAndTwistParameters);
		
		//use the updated chord and twist to generate a sequence of scaled and twisted aerofoils.
		var foils = this.updateBladeGeometry();	
		
		var bladeParameters = {
			surfaceCurves: foils,
			//vDegree: 3,
			vBspline: this.chordAndTwist.bspline.c, //chord and twist distributions as a bspline curve.
			wDegree: 3,
			alignedAtPercentChord: 0.,  //spanwise elements are aligned with leading edge.
			centroidAtPercentChord: undefined, //undefined: volume is centred at 'thickest' point on chord.		
		};
		
		this.initialise(bladeParameters); //initialise the NURBSSolid object.
		
		//shift the blade so that its long axis is along the midchord of the transition piece, not the leading edge.
		this.midhub_along_span_axis();	

		//Pitch the blade about its leading edge (currently Z-axis) by the pitch angle,
		//and rotate it so that the leading edge lies along the Y axis.
		//this.rotateAboutPoint(new Point_3D(0.,0.,0.), [-Math.PI/2., this.pitchAngle, -Math.PI/2.]); //[xy, xz, yz] radians.
		var pitchRadians = this.pitchAngle*Math.PI/180.;
		this.rotateAboutPoint(new Point_3D(0.,0.,0.), [-pitchRadians-(Math.PI/2.), 0., -Math.PI/2.]); //[xy, xz, yz] radians.	
		
		//experimental: raising degree recreates curves.
		//this.vDegree = 5;
		//this.degreeRaise_V();
	} 
	
	resetToOriginal(){
		var geometry = this.getNewcastleBladeDesign({numFoils:this.numFoils});	
		this.chordAndTwist = {};
		this.chordAndTwist.template = geometry.chordAndTwist;
		this.chordAndTwist.bspline = geometry.chordAndTwist_leastSquaresCurve;
	}
	
	getChordAndTwistDistribution(params){
		//input: (optional) params={chordAndTwist: LeastSquaresCurve() geometry definition object., isReset: bool}

		//ensure the this.chordAndTwist object exits.
		if(this.chordAndTwist === undefined)this.resetToOriginal();
		
		//default - used to retrieve current geometry (eg: for plotting).
		if(params === undefined)return this.chordAndTwist; 
		
		//reset
		if(params && params.isReset){
			this.resetToOriginal();
			return;
		}

		//does the blade exist, but without a properly defined chordAndTwist object?  
		//This occurs when older blades are loaded from file. 
		//Overwrite the default geometry definition with current one.
		if(params && params.isNew && this.curves){
			var chordAndTwist = this.getCurrentChordAndTwist(true); //get chord and twist distribution of surface curves.
			var chord = chordAndTwist.chordDistribution;
			var twist = chordAndTwist.twistDistribution;
			
			var c = new NURBSCurve();
			c.n = chord.length-1;
			c.p = this.q;
			c.U = [];
			for(var i=0;i<this.V.length;i++)c.U.push(this.V[i]);
			for(var i=0;i<=c.n;i++){
				c.P.push(new WeightedPoint_3D(chord[i].x, chord[i].y, twist[i].y, 1.));
			}

			this.chordAndTwist.bspline = new LeastSquaresCurve(); //dummy object.
			this.chordAndTwist.bspline.c = c;
		}
		
		return;
	}
	
	transitionHasChanged(){
		if(this.modelTransition !== this.hasTransition || 
		   this.transitionPercent !== this.currentTransitionPercent)return true;
		return false;
	}
	
	updateBladeGeometry(){
		var chordAndTwist = this.chordAndTwist.bspline.c;
		//chordAndTwist is a 3D B-spline curve with control points of: WeightedPoint_3D(r,chord,twist,weight)	
		
		//if the number of spanwise foils has changed, or if we're changing transition, revert to default geometry.
		if(this.numFoils != chordAndTwist.P.length || this.transitionHasChanged()){
			this.resetToOriginal();
			chordAndTwist = this.chordAndTwist.bspline.c;
		}
		
		//create an array of spanwise aerofoils based on these distributions.
		var foil = new Array(chordAndTwist.P.length); //SurfaceCurve[]
		
		var scalingFactor, vSpan, t;
		for(var i=0; i<chordAndTwist.P.length; i++){
			//make a copy of the template aerofoil
			if(i==0 && this.modelTransition)foil[i] = this.hubFoil.clone();
			else foil[i] = this.baseFoil.clone();

			//scale by chord length.
			scalingFactor = chordAndTwist.P[i].p.y*this.rBlade;
			//spanwise position.
			vSpan = chordAndTwist.P[i].p.x*this.rBlade;
			
			foil[i].initialise(scalingFactor, -chordAndTwist.P[i].p.z, vSpan, undefined, true);
		}
		
		//adjust the tip foil so that it extends right to the tip
		foil[foil.length-1].vSpan = this.rBlade;
		
		//update hubfoil so that its working curve includes the scaled/rotated 3D curve instead of 2D.
		if(this.modelTransition)this.hubFoil.c = foil[0].c.clone();	

		return foil;
	}
	
	midhub_along_span_axis(){
		//translate the blade so that the midchord of the connector piece is along the Y axis
		//instead of the leading edge
		var first = this.curves[0][0];
		first.findLeadingEdge();
		var v = first.LE.minus(first.c.P[0].p).multiply(0.5); //vector from LE to midchord.
		var j,k;
		for(var i=0;i<this.curves.length;i++){
			for(j=0;j<this.curves[i].length;j++){
				for(k=0;k<this.curves[i][j].c.P.length;k++){
					this.curves[i][j].c.P[k].p = this.curves[i][j].c.P[k].p.plus(v);
					this.curves[i][j].LE = this.curves[i][j].LE.plus(v);
				}
			}
		}
	}
	
	getCurrentChordAndTwist(justFoils){
		//return chord and twist distribution of surface aerofoil curves.
		
		var foils;
		if(justFoils){
			foils = this.curves[0]; //surface curves;
			for(var i=0;i<foils.length;i++){
				foils[i].findLeadingEdge();
				foils[i].TE = foils[i].c.P[0].p;
			}	
		}	
		else{
			foils = [];
			
			var numSlices = 100;
			var vStep = (this.V[this.V.length-1] - this.V[0])/(numSlices-1);

			var v, c;
			for(var i=0;i<numSlices;i++){
				v = i*vStep;
				
				c = this.extractSurfaceIsocurve(v);
				c.LE = SurfaceCurve.prototype.findLeadingEdge(c).point
				c.TE = c.P[0].p;
				
				foils.push(c);	
			}
		}
		
		var rBlade = this.curves[0][this.curves[0].length-1].c.P[0].p.y;
		
		var c, vec, rOnR, chordDistribution=[], twistDistribution=[];
		for(var i=0;i<foils.length;i++){
			rOnR = foils[i].TE.y/rBlade;
			
			vec = foils[i].TE.minus(foils[i].LE);

			chordDistribution.push({x:rOnR, y:vec.NORM()/rBlade});
			twistDistribution.push({x:rOnR, y:(-(Math.atan(vec.x/vec.z)*180./Math.PI)-this.pitchAngle)});
		}
		
		return {chordDistribution:chordDistribution, twistDistribution:twistDistribution};
	}
	
	getNewcastleBladeDesign(parameters){
		//input: parameters={numFoils}	
		// Chord distribution - 4th order poly fit to data in Anderson et al.
		// (1982) by DHW January 1997.
		
		//this routine is also called by the plotting routine to generate baseline geometry plot.
		if(!parameters || !parameters.numFoils){
			if(this.numFoils)parameters = {numFoils:this.numFoils};
			else{
				console.log("blade not properly defined.");
				return;
			}	
		}
		
		//discretisation of polynomial DHW design.
		var numPoints = parameters.numFoils*2;
		if(numPoints<40)numPoints=40;

		//create a data structure to store chord and twist distributions.  For convenience, chord and twist
		//are stored as arrays of WeightedPoint_2D, with .x being the span position (r) of the polynomial
		//distrubution, and y being normalised chord or twist.  The weight is used later for curve fitting.
		var design = {};
		design.chordAndTwist = new Array(numPoints); // each element is WeightedPoint_3D(r,chord,twist,weight);
		
		var rTip = 1.;//normalised blade length.
		var rHub = 0.06;//the start of the transition piece.  Rectangular connector piece to the left.
		//var rTransition = 0.288; //We want a maximum chord of 248mm on a 2.5 metre blade.  This occurs at r=~0.3.
		//if this.transitionPercent hasn't been defined, set it to a plausible default.  This occurs when old blade
		//designs are loaded from file and they don't have the this.transitionPercent parameter.
		this.transitionPercent = (this.transitionPercent !== undefined)?this.transitionPercent:28.8;
		var rTransition = this.transitionPercent/100.; //We want a maximum chord of 248mm on a 2.5 metre blade.  This occurs at r=~0.3.

		var rh = rTransition;
		var delr = (1. - rh)/(numPoints-2); //the first foil will be a transition piece.
		
		//1/12/2000 - new hub geometry.  The very first element will be the rectangular curve that
		//connects to the end piece.  So the first aerofoil section will be curve #2 (ie: index 1)
		
		//set chord and twist for first aerofoil.
		var first=0;
		if(this.modelTransition){
			//connector piece for 2.5 metre blade is 150mm per side.  So 150/2500 = 0.06.  Twist is 0 degrees.
			design.chordAndTwist[0] = new WeightedPoint_3D(rHub, 0.06, 0., 1.);
			
			first = 1;
		}

		this.hasTransition = this.modelTransition;
		this.currentTransitionPercent = this.transitionPercent;
		
		//the start of the blade proper.
		var r = rTransition, chord, twist;
		for(var i=first; i<numPoints; i++){
			//calculate chord at this span location
			chord =  2.5*(0.16165732 - 0.5847727*r+
							 1.0327255*Math.pow(r,2) -
							 0.8756711*Math.pow(r,3) +
							 0.2844545*Math.pow(r,4))/1.5;

			//calculate the twist at this span location
			if(r <= 0.7){
				twist = 54.16632 - r*(307.42939 - r*(719.549614 - r*(785.971096 - r*326.673372)));
			}
			else twist = 5.318999 - 7.059999*r;
			
			design.chordAndTwist[i] = new WeightedPoint_3D(r, chord, twist, 1.);		
			
			//next span position.
			r += delr;
		}

		/////
		///// fit a least-squares curve to this point data.
		/////

		//constrain end points so curve interpolates them.
		design.chordAndTwist[0].w = 0.;
		design.chordAndTwist[design.chordAndTwist.length-1].w = 0.;

		//derivative at start point.  Tack it onto the WeightedPoint_3D object - curve fitter can handle it.
	//	design.chordAndTwist[0].d = new WeightedPoint_3D(2, 0.45, 1., 0.);

		design.chordAndTwist[0].d = design.chordAndTwist[1].minus(design.chordAndTwist[0]);
		design.chordAndTwist[0].d.w = 5.;
		
	//	design.chordAndTwist[0].d = new WeightedPoint_3D(0.5, 0., 0., 0.);

		//fit a least-squares curve to this point data.
		var numFoils = parameters.numFoils, degree = 3; //defaults
		if(parameters.degree)degree=parameters.degree;
		
		design.chordAndTwist_leastSquaresCurve = new LeastSquaresCurve(design.chordAndTwist, parameters.numFoils-1, degree);

		return design;
	}
}
