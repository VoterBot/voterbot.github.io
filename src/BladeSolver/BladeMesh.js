"use strict";
class BladeMesh{
	constructor(x, meshSettings){
		//input (optional): BsplineBlade x, int b, int n, int m, int w, int t, double lambda, double u0, boolean s
		this.nPanels; //number of chordwise blade panels - int
		this.mPanels; //number of spanwise blade panels - int
		this.wPanels; //number of streamwise wake panels per wake revolution - int
		this.turns;    //number of revolutions of wake to model - int
		this.numBlades;//number of blades on rotor - int
		
		//some variables for the panel solver
		this.rBlade;	//float
		this.lambda;//the tip speed ratio - float
		this.omega;//the rotational speed of the rotor - float
		this.U_rotor;//wind speed at the rotor - float
		this.U_0;//free-stream wind speed - float
		this.isHelicalWake;//boolean
		this.allowWakeRollup;//boolean

		this.bladePanel; //an array of panels that define our blade mesh - SurfacePanel[]
		this.wakePanel;  //an array of wake panels - SurfacePanel[]

		this.element; // BladeElement[]
		
		if(x)this.setup(x, meshSettings);
	}	

	//methods
	setup(x, meshSettings){
		//input: BsplineBlade x, int b, int n, int m, int w, int t, double lambda, double u0, boolean s
		if(meshSettings){
			this.nPanels = meshSettings.nPanels;
			this.mPanels = meshSettings.mPanels;
			this.wPanels = meshSettings.wPanels;
			this.turns = meshSettings.turns;
			this.numBlades = meshSettings.numBlades;	
			//
			this.lambda = meshSettings.lambda;
			this.U_0 = meshSettings.U_0;
			this.isHelicalWake = meshSettings.isHelicalWake;
			this.allowWakeRollup = meshSettings.allowWakeRollup;
		}	

		this.rBlade = x.rBlade;      

		delete this.bladePanel; //attempt to force garbage collection.
		this.bladePanel = new Array(this.numBlades * this.nPanels * this.mPanels);
		for(var i=0;i<this.bladePanel.length;i++)this.bladePanel[i] = new SurfacePanel();
		
		delete this.wakePanel;
		this.wakePanel = new Array(this.numBlades * this.mPanels * this.wPanels * this.turns);
		for(var i=0;i<this.wakePanel.length;i++)this.wakePanel[i] = new SurfacePanel();
		
		delete this.element;
		this.element = new Array(this.numBlades * this.mPanels);
		for(var i=0;i<this.element.length;i++)this.element[i] = new BladeElement();	

		//create the mesh
		this.createMesh(x);
		//this.brickMesh = new Brick8Mesh(this);
	}  
	
	//copy constructor
	clone(a){
		//input: BladeMesh a (optional)
		//output: BladeMesh x
		
		if(a === undefined)a = this; //default copy constructor.
		
		var x = new BladeMesh();
		
		x.nPanels = a.nPanels;
		x.mPanels = a.mPanels;
		x.wPanels = a.wPanels;
		x.turns = a.turns;
		x.numBlades = a.numBlades;

		x.bladePanel = new Array(a.numBlades * a.nPanels * a.mPanels);
		for(var i=0; i<a.numBlades*a.nPanels*a.mPanels; i++)x.bladePanel[i] = SurfacePanel.prototype.clone(a.bladePanel[i]);

		x.wakePanel = new Array(a.numBlades * a.mPanels * a.wPanels * a.turns);
		for(var i=0; i<a.numBlades * a.mPanels * a.wPanels * a.turns; i++)x.wakePanel[i] = SurfacePanel.prototype.clone(a.wakePanel[i]);

		x.element = new Array(a.numBlades * a.mPanels);
		for(var i=0; i<a.numBlades * a.mPanels; i++)x.element[i] = BladeElement.prototype.clone(a.element[i]);

		x.rBlade = a.rBlade;
		x.lambda = a.lambda;
		x.omega = a.omega;
		x.U_0 = a.U_0;
		
		x.U_rotor = new Array(a.U_rotor.length);	
		for(var i=0;i<x.U_rotor.length;i++)x.U_rotor[i] = a.U_rotor[i];
		
		x.isHelicalWake = a.isHelicalWake;
		x.allowWakeRollup = a.allowWakeRollup;
		
		return x;
	}
	
	createMesh(x){
		//console.log("BladeMesh: createMesh()");
		//a routine to create a rectangular mesh from a b-spline blade representation	
		//input: BsplineBlade x
			
		if(!this.isHelicalWake)this.turns=1;
		
		//allocate space for our parameter arrays
		var u1 = [];//array of parameters in u direction (hubwise vertex strip on panel) - float[]
		var u2;//array of parameters in u direction (tipwise vertex strip on panel) - float[]
		var v = new Array(this.mPanels+1);//array of parameters in v direction	- float[]
				
		//**************
		//first create blade panels
		//**************	
		
		var foil = x.curves[0]; //array of surface curves.
		
		//find the U parameter of each blade foil leading edge
		for(var i=0;i<foil.length;i++)foil[i].findLeadingEdge();
		
		//use cosine-spacing to get mPanels span-wise panel parameters
		var theta, alpha, delta = Math.PI/this.mPanels;	//float.
		for(var i=0; i<=this.mPanels; i++){
			theta = i*delta;
			v[i] = (1. - Math.cos(theta))/2.;
		}
			
		//assemble panels by finding the 3D positions of their corner nodes
		var j, k = 0;
		for(var i=0; i<this.mPanels; i++){//do one span row of panels at a time
			//use cosine spacing to get nPanels chordwise
			//panel parameters for this foil strip
			u1 = this.get_U_parameters(x,v[i]);   //hubwise u parameters for this panel strip
			u2 = this.get_U_parameters(x,v[i+1]); //tipwise u parameters for this panel strip		
			
			//assemble foil and wake strip i
			for(j=0; j<this.nPanels; j++){
				if(j==0){
					this.bladePanel[k].node[3] = x.surfacePoint(u1[j], v[i]);
					this.bladePanel[k].node[2] = x.surfacePoint(u2[j], v[i+1]);
					//and set the TE wake nodes
					this.wakePanel[i].node[0] = this.bladePanel[k].node[3].clone();
					this.wakePanel[i].node[1] = this.bladePanel[k].node[2].clone();
					//tell this wake panel the indices of the TE blade panel that shed it
					this.wakePanel[i].upper = k; //upper surface of trailing edge
					this.wakePanel[i].lower = k+this.nPanels-1;	
				}
				else{
					this.bladePanel[k].node[3] = this.bladePanel[k-1].node[0].clone();
					this.bladePanel[k].node[2] = this.bladePanel[k-1].node[1].clone();
				}
				if(j<this.nPanels-1){
					this.bladePanel[k].node[1] = x.surfacePoint(u2[j+1], v[i+1]);
					this.bladePanel[k].node[0] = x.surfacePoint(u1[j+1], v[i]);
				}
				else{//Lower TE panel
					this.bladePanel[k].node[1] = this.bladePanel[k-this.nPanels+1].node[2].clone();
					this.bladePanel[k].node[0] = this.bladePanel[k-this.nPanels+1].node[3].clone();
					if(this.bladePanel[k].node[0].u<0.5)this.bladePanel[k].node[0].u = 1.-this.bladePanel[k].node[0].u;
					if(this.bladePanel[k].node[1].u<0.5)this.bladePanel[k].node[1].u = 1.-this.bladePanel[k].node[1].u;				
				}
				this.bladePanel[k].m_panel = i; //tag the panel with its element number
				this.bladePanel[k].theta = 0.0; //Key blade is aligned with global Y axis
				this.bladePanel[k].r_on_R = (this.bladePanel[k].node[3].y +
				   0.5*(this.bladePanel[k].node[2].y - this.bladePanel[k].node[3].y))/x.rBlade;
				k++;
			}
		}
		
		//if the rotor has more than one blade, make (numBlades-1) additional copies of the above blade mesh
		if(this.numBlades > 1)this.copyBladeMesh();
		
		//**************
		//create a default wake geometry - using eqs 40 from Preuss, Suciu and Morino (1980)
		//**************	
		if(this.lambda < 1.E-6)this.lambda=1.E-6;
		this.omega  = this.lambda*this.U_0/x.rBlade;
		this.U_rotor = new Array(this.mPanels);
		for(var i=0;i<this.mPanels;i++)this.U_rotor[i] = this.U_0;
		this.makeWake();
		
		//new: create a mesh for the FE solver using 8-node bricks
		//create_FEA_mesh(mesh);
		
		//if the potential flow solver has been turned off by the user, the panel outward
		//normals and local coordinate systems won't be calculated, so won't be available
		//to view on our nice graphical display.  Let's calculate them here.
		//
		//blade panels
		for(var i=0;i<this.bladePanel.length;i++){
			this.bladePanel[i].getHyperboloidalCornerPoints();
			this.bladePanel[i].calculatePanelOutwardNormal();
		}
		//wake panels
		for(var i=0;i<this.wakePanel.length;i++){
			this.wakePanel[i].getHyperboloidalCornerPoints();
			this.wakePanel[i].calculatePanelOutwardNormal();
		}

		return;
	}
	
	get_U_parameters(x, v){
		//find nPanels+1 U parameters using cosine spacing
		//input: BsplineBlade x, double v	
		//input:  x (Bspline blade representation)
		//input:  v (spanwise parameter for this row of panel vertices)
		//output: u[] (array of u parameters for this strip of panel vertices)

		var foil = x.curves[0]; //array of surface curves.
		
		//first, find the two foils surrounding this v parameter
		var j=0;
		for(j=0;j<foil.length;j++){
			if((foil[j].vSpan-foil[0].vSpan)/(x.rBlade-foil[0].vSpan) > v)break;
		}
		
		if(j>foil.length-1)j=foil.length-1;

		//now, interpolate linearly between the U_LE paramters of these two
		//foils to get the U_LE parameter for this vertex strip
		var a = (foil[j-1].vSpan-foil[0].vSpan)/(x.rBlade-foil[0].vSpan); //float
		var b = (foil[j].vSpan-foil[0].vSpan)/(x.rBlade-foil[0].vSpan); //float
		
		var ratio = (v-a)/(b-a); //float
		
		//find the linear interpolation LE U parameter for this vertex strip
		var u_LE = foil[j-1].u_LE + ratio*(foil[j].u_LE - foil[j-1].u_LE); //float
		
		//use cosine-spacing to get nPanels chord-wise panel parameters
		var u = new Array(this.nPanels+1); //float[]
		var theta, delta = 2.*Math.PI/this.nPanels;
		for(var i=0; i<=this.nPanels; i++){
			theta = i*delta;

			//half of the panels go on upper half of blade U={0,U_LE}
			if(theta <= Math.PI)u[i] = u_LE*(1.-Math.cos(theta))/2.;
			//and half on the bottom U=[U_LE,1}
			else u[i] = u_LE + (1.-u_LE)*(1.-Math.cos(theta - Math.PI))/2.;
		}

		return u;
	}
	
	copyBladeMesh(){
		//make a copy of the blade mesh for each additional blade on the rotor
		var numPanels = this.nPanels * this.mPanels;

		//determine the angle by which the new blade(s) will be offset (about rotor axis)
		var theta = 2.*Math.PI/this.numBlades;

		//copy the new blade mesh panels into their new positions
		var costh, sinth, y, z, i, j;
		for(var blade = 1; blade<this.numBlades; blade++){
			costh = Math.cos(blade*theta);
			sinth = Math.sin(blade*theta);
			for(i=0; i<numPanels; i++){
				for(j=0;j<4;j++){
					y = this.bladePanel[i].node[j].y;
					z = this.bladePanel[i].node[j].z;
					//rotate about the global x axis
					this.bladePanel[i+(blade*numPanels)].node[j].x = this.bladePanel[i].node[j].x;
					this.bladePanel[i+(blade*numPanels)].node[j].y = y*costh - z*sinth;
					this.bladePanel[i+(blade*numPanels)].node[j].z = y*sinth + z*costh;
				}
				this.bladePanel[i+(blade*numPanels)].m_panel = this.bladePanel[i].m_panel;
				this.bladePanel[i+(blade*numPanels)].theta   = blade*theta;
				this.bladePanel[i+(blade*numPanels)].r_on_R  = this.bladePanel[i].r_on_R;
			}
		}

		return;
	}
	
	makeWake(){	
		//**************
		//create wake panels - using eqs 40 from Preuss, Suciu and Morino (1980)
		//**************	
		
		var numStreamwisePanels = this.turns*this.wPanels;//helical wake
		
		//if we're not allowing the wake to rollup due to the distribution of U_rotor along the blade, 
		//use an average value of U_rotor.
		var U_rotor;
		if(!this.allowWakeRollup)U_rotor = this.get_mean_U_rotor();
		
		var k = 0, i, j;
		var angle=0., alpha=0.;
		for(var i=0; i<numStreamwisePanels; i++){ //do one wake strip at a time (parallel to TE)
			//set initial wake step
			alpha=0.;
			if(this.isHelicalWake){
			   //16/5/2001 - shorten the wake panels closest to the TE (ie: the first turn)
			   if(i<this.wPanels)alpha = 2.0*Math.PI/this.wPanels/this.wPanels;//shortened helical wake
			   else alpha = 2.0*Math.PI/this.wPanels;//helical wake		
			}
			else  alpha = 1.0; //straight wake
			
			for(j=0; j<this.mPanels; j++){
				if(this.allowWakeRollup)U_rotor = this.U_rotor[j];
				if(i==0)this.wakePanel[k] = this.getFirstWakePanel(U_rotor, this.wakePanel[k],alpha);
				else{
					this.wakePanel[k].node[0] = this.wakePanel[k-this.mPanels].node[3].clone();
					this.wakePanel[k].node[1] = this.wakePanel[k-this.mPanels].node[2].clone();
					//and tell this wake panel which pair of TE blade panels shed it
					this.wakePanel[k].upper = this.wakePanel[k-this.mPanels].upper;
					this.wakePanel[k].lower = this.wakePanel[k-this.mPanels].lower;
			
					this.wakePanel[k].node[2] = this.wakePoint(U_rotor, this.wakePanel[j].node[2], angle);
					
					/////Modification to simple wake geometry:  instead of having independent streamwise
					//wake strips, tack node[3] of each spanwise panel onto node[2] of the panel beside it.
					if(j==0)this.wakePanel[k].node[3] = this.wakePoint(U_rotor, this.wakePanel[j].node[3], angle);
					else this.wakePanel[k].node[3] = this.wakePanel[k-1].node[2].clone();
				}

				k++;
			}
			angle+=alpha; //incrementing the angle here produces one extra near-wake panel		
		}
		
		//make a copy of the above wake geometry for each blade on a multi-bladed rotor
		if(this.numBlades > 1)this.copyWakeMesh();

		return;
	}
	
	getFirstWakePanel(U_rotor, x, alpha){
		//input: float U_rotor, SurfacePanel x, float alpha
		//output: SurfacePanel
		
		var panel = x.clone();

		//first, find the point that lies half-way between the upper and lower
		//collocation points on the TE blade panels that shed this wake panel
		var A = this.bladePanel[panel.upper].node[0].plus(
					this.bladePanel[panel.lower].node[3].minus(
						this.bladePanel[panel.upper].node[0]).multiply(0.5)); //Point_3D

		var B = this.bladePanel[panel.upper].node[3]; //Point_3D
		var BminusA = B.minus(A);//Point_3D
		var direction  = BminusA.divide(BminusA.NORM());//Point_3D
		//29/3/2001: this direction causes problems near the hub where the sections are not
		//aerofoils - the first wake panel here makes almost a 90 degree angle with the TE!
		//this is due to the rounded edge - the routine is working like it should.
		panel.node[3] = this.firstWakePoint(U_rotor, panel.node[0], alpha, direction);

		A = this.bladePanel[panel.upper].node[1].plus(
				(this.bladePanel[panel.lower].node[2].minus(
					this.bladePanel[panel.upper].node[1])).multiply(0.5));

		B = this.bladePanel[panel.upper].node[2];
		BminusA = B.minus(A);		
		direction = BminusA.divide(BminusA.NORM());

		panel.node[2] = this.firstWakePoint(U_rotor, panel.node[1], alpha, direction);

		return panel;		
	}
	
	firstWakePoint(U_rotor, TE, alpha, direction){
		//input: Point_3D TE, double alpha, Point_3D direction
		//output: Point_3D

		var alpha_TE = Math.atan(TE.z/TE.y);
		var R_TE = Math.sqrt(TE.y*TE.y + TE.z*TE.z);
		
		var point = new Point_3D();
		if(this.isHelicalWake){ //helical wake
			//eqs 40 in Preuss et. al. 1980
			//point.x = TE.x + U_0*alpha/pow(OMEGA,1.225);
			point.x = TE.x + U_rotor*alpha/this.omega;
			//point.x = TE.x + alpha*U_0/OMEGA;	
			point.y = R_TE*Math.cos(alpha_TE - alpha);
			point.z = R_TE*Math.sin(alpha_TE - alpha);

			//adjust the position of the first wake points so that
			//wake is shed smoothly from the TE
			//20/7/2000
			//point = TE + 0.5*direction*sqrt(TE.x*TE.x + TE.z*TE.z);
			//point = TE + direction*0.01;
			//pick the shortest of the two first panel choices
			if(point.minus(TE).NORM() > 0.0375*Math.sqrt(TE.x*TE.x + TE.z*TE.z))
				point = TE.plus(direction.multiply(0.0375*Math.sqrt(TE.x*TE.x + TE.z*TE.z)));
			else	
				point = TE.plus(direction.multiply(point.minus(TE).NORM()));
		}
		else{ //straight wake for testing
		
	//		   point.x = TE.x + alpha*WAKE_LENGTH*hub_chord/wPanels;
	//		   point.y = TE.y;
	//		   point.z = TE.z;

		   //point = TE + 0.1*direction*NORM(TE-point);
		}

		return point;
	}
	
	wakePoint(U_rotor, TE, alpha){
		//input: Point_3D TE, float alpha
		
		var alpha_TE = Math.atan(TE.z/TE.y);
		var R_TE = Math.sqrt(TE.y*TE.y + TE.z*TE.z);
		
		var point = new Point_3D();
		if(this.isHelicalWake){ //helical wake
			//eqs 40 in Preuss et. al. 1980


			//find the down-stream position of this wake point.  If you let the
			//wake travel at the speed of the undisturbed flow (U0), you get:	
			//point.x = TE.x + U_0*alpha/OMEGA;
			//but this doesn't give very good Cp results, mostly because it
			//doesn't take the blade thrust into account - ie: the higher the
			//thrust, the slower the downstream flow will be, so the wake
			//will be closer to the rotor, and thus influence it more.  The
			//most accurate way to find the speed of the downstream flow is to
			//take the rotor thrust into consideration, which means that
			//you have to iterate - you have to solve the linear system
			//as many times as the iteration takes to converge.  This is
			//a pain in the arse - it takes far too long.

			//or how about basing the pitch on the flow through the rotor?  
			point.x = TE.x + U_rotor*alpha/this.omega;

			point.y = R_TE*Math.cos(alpha_TE - alpha);
			point.z = R_TE*Math.sin(alpha_TE - alpha);

			//account for wake expansion
			//delta_r = sqrt((point.x/R_blade)*LAMBDA/100.0);
			//point.y = point.y*(1.+delta_r);
			//point.z = point.z*(1.+delta_r);
		}
		else{ //straight wake for testing
		
	//		   point.x = TE.x + alpha*WAKE_LENGTH*hub_chord/wPanels;
	//		   point.y = TE.y;
	//		   point.z = TE.z;
		   
		}

		return point;
	}
	
	get_mean_U_rotor(){
		//output: float

		//find the average value of U_rotor, weighted by the length of each trailing edge panel.
		var upper, lower;
		var length, total_length = 0., mean_U_rotor = 0.;
		for(var i=0; i<this.mPanels; i++){
			upper = this.wakePanel[i].upper;
			lower = this.wakePanel[i].lower;

			//ignore the hubwise 20% of the blade: 
			if(this.bladePanel[upper].r_on_R > 0.2){
				length = this.bladePanel[lower].node[3].minus(this.bladePanel[lower].node[2]).NORM();	   
				mean_U_rotor += this.U_rotor[i]*length;
				total_length += length;
			}	
		}

		return (mean_U_rotor/total_length);		
	}
	
	copyWakeMesh(){
		//make a copy of the blade mesh for each additional blade on the rotor

		var numBladePanels = this.mPanels * this.nPanels;	
		var numWakePanels = this.mPanels * this.wPanels * this.turns;

		//determine the angle by which the new wake panel will be offset (about rotor axis)
		var theta = 2.*Math.PI/this.numBlades;

		//copy the new blade mesh panels into their new positions
		var costh, sinth, i, j, y, z;
		for(var blade = 1; blade<this.numBlades; blade++){
		  costh = Math.cos(blade*theta);
		  sinth = Math.sin(blade*theta);
		  for(i=0; i<numWakePanels; i++){
			 for(j=0;j<4;j++){
				y = this.wakePanel[i].node[j].y;
				z = this.wakePanel[i].node[j].z;
				this.wakePanel[i+(blade*numWakePanels)].node[j].x = this.wakePanel[i].node[j].x;
				this.wakePanel[i+(blade*numWakePanels)].node[j].y = y*costh - z*sinth;
				this.wakePanel[i+(blade*numWakePanels)].node[j].z = y*sinth + z*costh;
			 }
			 //determine the TE blade elements that emitted this wake element
			 this.wakePanel[i+(blade*numWakePanels)].upper = this.wakePanel[i].upper + numBladePanels;
			 this.wakePanel[i+(blade*numWakePanels)].lower = this.wakePanel[i].lower + numBladePanels;
		   }
		}

		return;
	}
	
	outputCornerPoints(){
		var B = this.nPanels*this.mPanels;
		
		var outfile="";

		outfile+=("\n------------------\n");
		outfile+=("nPanels: "+this.nPanels+" mPanels: "+this.mPanels+"\n");
		outfile+=("Total Panels: "+B);
		outfile+=("\n------------------\n");
		
		var p0, p1, p2, p3;
		for(var i=0; i<B; i++){
			p0 = this.bladePanel[i].node[0];
			p1 = this.bladePanel[i].node[1];
			p2 = this.bladePanel[i].node[2];
			p3 = this.bladePanel[i].node[3];
			outfile+=("panel["+i+"].node: [0]("+p0.x+","+p0.y+","+p0.z+") [1]("+p1.x+","+p1.y+","+p1.z+") [2]("+p2.x+","+p2.y+","+p2.z+") [3]("+p3.x+","+p3.y+","+p3.z+")\n");
		}		
				
		outfile+=("\n------------------\n");
		
		//saveStrings("mesh_bladePanels.dat",outfile);
		console.log(outfile);
	}
	
	outputWakeCornerPoints(){
		var W = this.numBlades * this.mPanels * this.wPanels * this.turns;
			
		var outfile="";

		outfile+=("\n------------------\n");
		outfile+=("numBlades="+this.numBlades+" mPanels: "+this.mPanels+" wPanels: "+this.wPanels+" turns: "+this.turns+"\n");
		outfile+=("Total Panels: "+W);
		outfile+=("\n------------------\n");
		
		var p0, p1, p2, p3;
		for(var i=0; i<W; i++){
			p0 = this.wakePanel[i].node[0];
			p1 = this.wakePanel[i].node[1];
			p2 = this.wakePanel[i].node[2];
			p3 = this.wakePanel[i].node[3];
			outfile+=("panel["+i+"].node: [0]("+p0.x+","+p0.y+","+p0.z+") [1]("+p1.x+","+p1.y+","+p1.z+") [2]("+p2.x+","+p2.y+","+p2.z+") [3]("+p3.x+","+p3.y+","+p3.z+")\n");
		}		
				
		outfile+=("\n------------------\n");
		
		//saveStrings("mesh_wakePanels.dat",outfile);
		console.log(outfile);
	}
	
	solve(resultsHandler){
		//input: (optional) object with functions to handle results.
		//run the panel method solver on this mesh.
		
		if(resultsHandler === undefined)resultsHandler = {};
		
		//////////
		//	The stand-alone function.
		//	var solver = new PotentialFlow(this);
		//////////

		///////
		// The web worker.
		var worker = new Worker('PotentialFlow.js');
		worker.text = {
			permanent: "",
			temp: "",
			result: "",
		};
		
		worker.onmessage = function(e){
			if(e.data.solvedMesh){
				//the object returned is a JSON-stringified version of a BladeMesh.  Return the full object.		
				var x = BladeMesh.prototype.clone(e.data.solvedMesh);
				if(resultsHandler.updateMeshDisplay)resultsHandler.updateMeshDisplay(x);
			}
			if(e.data.message){
				if(e.data.isWipe)this.text.temp = "";		
				if(typeof e.data.message === "string"){
					if(e.data.isPermanent)this.text.permanent += e.data.message+"<br>";					
					else this.text.temp += e.data.message+"<br>";
				}
				else if(e.data.message.length > 0){
					this.text.result = "";
					for(var i=0;i<e.data.message.length;i++)this.text.result += e.data.message[i]+"<br>";
				}
					
				if(resultsHandler.setSolverResults){
					resultsHandler.setSolverResults(this.text.permanent + this.text.result + this.text.temp);
				}	
				console.log(e.data.message);
			}
		};
		worker.onerror = function(e) {
			alert('Error: Line ' + e.lineno + ' in ' + e.filename + ': ' + e.message);
		};

		//start the worker
		worker.postMessage({'cmd': 'PotentialFlow', 'bladeMesh': this,});
		
		return worker;
	}
	
	stopSolver(worker){
		if(worker){
			console.log("STOPPING SOLVER ...");
			worker.terminate();
			worker = undefined;
			console.log("STOPPED.");	
		}	
		
		return worker;
	}
}
//////
//////
class SurfacePanel{
	constructor(){
		this.node;	//the corner nodes of this quadrilateral panel - Point_3D[]
		this.P;		//corner points in hyperboloidal form - Point_3D[] 
					//P[0] is the panel centre point.  Used as collocation point.
		this.n;		//the unit outward normal of this panel - Point_3D 
		this.V_w;	//unperturbed wind velocity in rotating frame at this panel - Point_3D 
		this.dphi_dn;//the zero-flowthrough velocity boundary condition - float
		this.phi;	//the value of the velocity potential at the centre of this panel - float
		this.upper;  //for wake panels, the indices of the upper and lower - int
		this.lower;	//blade trailing-edge panels from which this panel emanated - int 
		this.m_panel;//span-wise blade element to which this panel belongs - int
		this.a1;		//Point_3D 
		this.a2;		//covariant base vectors along panel's local axes (xi and eta) - Point_3D 
		this.unit_Xi;//unit-vectors defining the local coordinate system - Point_3D 
		this.unit_Eta;//Point_3D 
		//Point_3D c_a1,c_a2,c_n;	//contravariant base vectors
		this.x_on_C;		//the chord-span position of this panel's centre point - float
		this.Cp;		//the pressure coefficient on this panel - float
		this.pressure;// float
		this.force;	// float
		this.theta;	//the angle that the entire blade that this panel belongs to
					//makes in the global Y-Z plane (measured counterclockwise from
					//the positive Y axis ie: positive RH rotation about X axis) - float
		this.r_on_R;		//spanwise position of this panel's centroid - float
		
		this.initialise();
	}	

	//methods
	initialise(){
		this.node = new Array(4);
		for(var i=0;i<this.node.length;i++)this.node[i] = new Point_3D();

		this.P = new Array(4);
		for(var i=0;i<this.P.length;i++)this.P[i] = new Point_3D();
		
		this.n = new Point_3D();
		this.V_w = new Point_3D();
		this.dphi_dn = 0.;
		this.phi = 0.;
		this.upper = 0;
		this.lower = 0;
		this.m_panel = 0;
		this.a1 = new Point_3D();
		this.a2 = new Point_3D();		
		this.unit_Xi = new Point_3D();
		this.unit_Eta = new Point_3D();
		this.x_on_C = 0.;
		this.Cp = 0.;
		this.pressure = 0.;
		this.force = 0.;
		this.theta = 0.;
		this.r_on_R = 0.;
	}	
	
	clone(a){
		//input: SurfacePanel a (optional)
		//output: SurfacePanel x
		
		if(a === undefined)a = this; //default copy constructor.
		
		var x = new SurfacePanel();

		x.node = new Array(a.node.length);
		for(var i=0;i<a.node.length;i++)x.node[i] = Point_3D.prototype.clone(a.node[i]);

		x.P = new Array(a.P.length);
		for(var i=0;i<a.P.length;i++)x.P[i] = Point_3D.prototype.clone(a.P[i]);
		
		x.n = Point_3D.prototype.clone(a.n);
		x.V_w = Point_3D.prototype.clone(a.V_w);
		x.dphi_dn = a.dphi_dn;
		x.phi = a.phi;
		x.upper = a.upper;
		x.lower = a.lower;
		x.m_panel = a.m_panel;
		x.a1 = Point_3D.prototype.clone(a.a1);
		x.a2 = Point_3D.prototype.clone(a.a2);		
		x.unit_Xi = Point_3D.prototype.clone(a.unit_Xi);
		x.unit_Eta = Point_3D.prototype.clone(a.unit_Eta);
		x.x_on_C = a.x_on_C;
		x.Cp = a.Cp;
		x.pressure = a.pressure;
		x.force = a.force;
		x.theta = a.theta;
		x.r_on_R = a.r_on_R;
		
		return x;
	}
	
	getHyperboloidalCornerPoints(){
		//equation 12 in Morino et. al. (1975)
		
		//number nodes CLOCKWISE to give a RH coordinate system
		this.P[0] = (this.node[2].plus(this.node[3]).plus(
					this.node[1]).plus(this.node[0])).multiply(0.25);
		this.P[1] = (this.node[2].plus(this.node[3]).minus(
					this.node[1]).minus(this.node[0])).multiply(0.25);
		this.P[2] = (this.node[2].minus(this.node[3]).plus(
					this.node[1]).minus(this.node[0])).multiply(0.25);
		this.P[3] = (this.node[2].minus(this.node[3]).minus(
					this.node[1]).plus(this.node[0])).multiply(0.25);

		return;
	}
	
	calculatePanelOutwardNormal(){
		//first, find the base vectors for this panel
		this.get_a1(0);		
		this.get_a2(0);

		//calculate the panel outward unit normal from a1 and a2
		this.get_n(this.a1, this.a2);

		//calculate the unit_Xi and unit_Eta vectors for this panel
		this.unit_Xi = this.a1.divide(this.a1.NORM());
		this.unit_Eta = this.a2.divide(this.a2.NORM());

		return;
	}
	
	get_a1(eta){
		//eq 18 Morino
		//input: float eta
		this.a1 = this.P[1].plus(this.P[3].multiply(eta));	
	}
	get_a2(xi){
		//eq 19 Morino
		//input: float xi
		this.a2 = this.P[2].plus(this.P[3].multiply(xi));
	}	
	get_n(a1, a2){
		//eq 20 Morino
		//input: Point_3D a1, Point_3D a2
		var B = a1.cross(a2); //Point_3D
		this.n = B.divide(B.NORM());
	}
}	
/////
/////
class BladeElement{
	constructor(){
		this.chord;	//blade element chord
		this.twist;   //blade element twist (degrees)
		this.alpha;	//blade element angle of attack (degrees)
		this.theta;	//blade element angle of incidence (degrees)
		this.Re;	//local Reynolds number at LE of element
		this.lift;
		this.drag;
		this.Cl;
		this.Cd;
	}	
	
	//methods
	clone(a){
		//input: BladeElement a (optional)
		//output: BladeElement x
		
		if(a === undefined)a = this; //default copy constructor.
		
		var x = new BladeElement();
		
		x.chord=a.chord;
		x.twist=a.twist;
		x.alpha=a.alpha;
		x.theta=a.theta;
		x.Re=a.Re;
		x.lift=a.lift;
		x.drag=a.drag;
		x.Cl=a.Cl;
		x.Cd=a.Cd;
		
		return x;
	}
}	
////
/*
function testBladeMesh(){
	var bladeSettings = {
		numFoils: 15,
		rBlade: 2.5,	
		pitchAngle: 0.,
		modelTransition: true, //model the transition section: boolean.
	};	
	
	var baseAerofoil = new SurfaceCurve();
	baseAerofoil.initialiseFromFile(loadData("7062_bspline"));
	var hubAerofoil = new SurfaceCurve();
	hubAerofoil.initialiseFromFile(loadData("rectangle_bspline"));
	
	var x = new BsplineBlade();
	x.createBlade(bladeSettings, baseAerofoil, hubAerofoil);

	var meshSettings = {
		numBlades: 2,
		nPanels: 30,
		mPanels: 10,
		wPanels: 20,
		turns: 10,
		lambda: 8.5,
		U_0: 10,
		isHelicalWake: true,
	}

	var mesh = new BladeMesh(x, meshSettings);
	mesh.solve();

	//mesh.outputCornerPoints();
	//mesh.outputWakeCornerPoints();		
}

function loadData(fileID){
	//extract data from the datafile loaded into an iframe in the HTML.
	//NOTE:  To run locally, Chrome requires this switch on the command line: 
	//  --allow-file-access-from-files
	//Otherwise, it gets angry at you trying to load data from a local iframe.
	//Should work fine when served by a web server.
	var oFrame = document.getElementById(fileID);
	return oFrame.contentWindow.document.body.childNodes[0].innerHTML;
};
*/