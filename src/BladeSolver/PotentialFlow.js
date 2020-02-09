"use strict";
class PotentialFlow{
	//Numerical panel method solver.
	//Routine can be invoked as a function or as a message-passing web worker.
	constructor(bladeMesh, isWorker){	
		//some physical constants
		this.rho = 1.204; //density of air at 20 degrees C (kg/m^3)
		this.visc = 1.151163E-5;//kinematic viscosity of air at 20 degrees C (m^2/s)
		this.rho_gfrp = 1500.0; //kg per cubic metre
		this.rho_foam = 150.0;  //kg per cubic metre
		
		//
		this.power;
		this.Cp; //power coefficient
		this.thrust;
		this.Cf; //thrust coefficient
		this.torque;
		
		//has this routine been invoked as a function or a stand-alone web worker?
		this.isWorker = isWorker;
		if(isWorker)self.importScripts('Point_3D.js','D_matrix.js','BladeMesh.js');
		
		//solve.
		var x = this.solve(bladeMesh);
		//return solution.
		if(x && this.isWorker)self.postMessage({'hasSolved': true, 'solvedMesh': x});
		else return x;
	}

	//methods
	solve(bladeMesh){
		//Solve Laplace's equation at the collocation point of each panel to get potential distribution.
		//Use Bernoulli's theorem to calculate pressure distribution, then torque, thrust and power.
		
		//If this is a web worker the bladeMesh object passed to it will be a
		//JSON-stringified object with no methods.  Recreate the full BladeMesh object.
		//This works just as well when invoked as a function.
		var x = BladeMesh.prototype.clone(bladeMesh);

		var B = x.nPanels * x.mPanels;
		var W = x.mPanels * x.turns * x.wPanels;
		
		var i, j;

		//output number of blade panels, and tag the text as permanent.
		this.resultsOutput("num_blade_panels: "+(x.numBlades*B), null, true);
		this.resultsOutput("num_wake_panels: "+(x.numBlades*W), null, true);
		
		//find the corner nodes of each blade panel in hyperboloidal form
		//(eq (12) in Morino et. al.1975)
		for(i=0;i<x.bladePanel.length;i++){
			x.bladePanel[i].getHyperboloidalCornerPoints();
			x.bladePanel[i].calculatePanelOutwardNormal();
		}
		
		//find blade panel surface area for future force/power calculations.
		this.calculateBladePanelSurfaceArea(x);

		//form influence matrices (blade panels)
		var Chk = this.get_Chk_and_Bhk(x, true);
		var Bhk = this.get_Chk_and_Bhk(x, false);
		
		//We must create new wake geometry (and new wake influence coefficients) for
		//the LAMBDA we are interested in studying
		//helical wake - turbine blade
		var lambda = x.lambda;
		var omega = x.lambda*x.U_0/x.rBlade;
		var nBladeIterations = 4; 
		var nWakeIterations = 15; //number of iterations for rotor geometry		   
		if(!x.isHelicalWake){
			lambda = 0.0;
			omega  = 0.0;
			nBladeIterations = 1; 
			nWakeIterations = 1; //number of iterations for straight geometry
		}

		this.Cp=0.;//power coefficient
		x.U_rotor = new Array(x.mPanels); //U_rotor at each spanwise position.
		var U_rotor_old={U_minus: new Array(x.mPanels), U_minus_minus: new Array(x.mPanels)};		
		for(i=0;i<x.mPanels;i++){
			x.U_rotor[i] = x.U_0; //a first guess at U_1 for the first iteration.
			U_rotor_old.U_minus[i] = x.U_0; //a first guess at U_1 for the first iteration.			
		}	
		//create a new wake for this LAMBDA and the current U_rotor
		var h, sigma, RHS, Thm, Chk_a, phi, blade, h, dP, wakeConvergence, iteration;		
		for(var wakeIteration=1;wakeIteration<=nWakeIterations;wakeIteration++){
			this.resultsOutput(("wake iteration "+wakeIteration+" of "+nWakeIterations), true);
		
			x.makeWake();
	
			//calculate the hyperboloidal corner points of wake panels
			//Morino et. al. (1975) (eq 12)
			for(i=0;i<x.wakePanel.length;i++){
				x.wakePanel[i].getHyperboloidalCornerPoints();
				x.wakePanel[i].calculatePanelOutwardNormal();
			}
			
			//calculate the boundary conditions -> source strengths
			x = this.boundary_conditions(x);
		
			//correct the boundary conditions for viscous effects
			x = this.get_elemental_quantities(x);
			x = this.viscous_boundary_correction(x);
	
			//load the source-strengths into the dphi_dn column vector
			sigma = new D_matrix(B,1);//a column vector containing the known panel BC source value		    
			for(h=0; h<B; h++)sigma.A[h][0] = x.bladePanel[h].dphi_dn;

			//multiply the Bhk matrix with the dphi_dn vector
			RHS = Bhk.cross(sigma).negative();  //-(B x sigma) provides RHS of eq 12, Preuss et.al.1980

			//account for wake panel influences on blade
			//apply_Fhm(X, Chk_a);
			Thm = this.get_Thm(x);

			for(iteration=1; iteration<=nBladeIterations; iteration++){
				this.resultsOutput(("pressure equalisation iteration "+iteration+" of "+nBladeIterations), true);
				
				//apply wake influence to blade influence matrix
				Chk_a = Chk.plus(Thm);
				
				//solve eq 12 in Preuss et. al. 1980 to get the unknown
				//velocity potential on each panel
				phi = Chk_a.solve(RHS); //solve the linear system (see D_matrix for details)
	
				//copy these values of phi back to the collocation point of each panel
				for(blade=0;blade<x.numBlades;blade++)
				   for(h=0; h<B; h++)x.bladePanel[h+(blade*B)].phi = phi.A[h][0];	    	
				
				//simple pressure determination - see Katz and Plotkin
				x = this.get_pressure(x);
	
				//adjust the wake strength to achieve zero-pressure difference at TE  	
				dP = this.get_dP(x);
				Thm = this.iterate_Thm(x, Thm, dP, iteration);
			}//end of Kutta condition iteration loop
			
			//output_delta_phi(x);//output the circulation on the blade (for plotting)	
		
			//calculate the power, thrust, power coefficient etc. for this blade
			this.calculatePower(x);
			
			//wake model adjustment: use circulation around the blade to find wind speed at rotor.
			x = this.find_U_rotor(x);

			//finally, calculate lift and drag for each blade element.
			x = this.calculateLiftAndDrag(x);
			////
			////
	
			//check for wake pitch convergence
			wakeConvergence = this.wake_convergence(x, U_rotor_old, wakeIteration);
			//output results of this iteration.
			this.outputResults(wakeIteration, nWakeIterations, wakeConvergence, x);

			//stop iteration if wake has converged.
			if(wakeConvergence.isConverged){
				break; //wake pitch isn't changing very fast - converged!
			}
			//make a copy of previous U_rotor iterations.
			for(i=0;i<x.mPanels;i++){
				U_rotor_old.U_minus_minus[i] = U_rotor_old.U_minus[i];			
				U_rotor_old.U_minus[i] = x.U_rotor[i];
			}	
		}//end of wake_iteration loop

		//approximate solid blade geometry with eight-node bricks for FE analysis
		//create_FEA_mesh(X);

		/*
		//23/6/2001: output panel centre coords with pressure data for Phil
		output_pressure_plot(X, Cp);
		
		output_U1_vs_tsr(U_rotor,(int)tsr);
	
		//output torque and force along blade	
		output_torque_vs_span(X);
		
		//output lift coefficient along blade
		output_lift_vs_span(X);
		*/		
		
		/*
			resultsOutput.setText("running FEA solver.  Please wait....\n");
			//run the slffea program
			time(&tv1);	
			system("/home/markh/downloads/FE_stuff/slffea-1.1/brick/brick/br fea_output/blade");
			time(&tv2);
			resultsOutput.setText("total FEA code run time = %d seconds\n",tv2 - tv1);	
		*/	
		return x;
	}
	
	outputResults(wakeIteration, nWakeIterations, wakeConvergence, x){
		var output = [];

		output.push("-----------<br>wake iteration: "+wakeIteration+" of "+nWakeIterations+"<br>-----------");		
		output.push(wakeConvergence.oldWakeResults);
		output.push("Cp: "+this.Cp.round(4));
		output.push("U_0: "+x.U_0.round(1)+" m/s");
		//output.push("U_rotor: "+x.U_rotor);
		output.push("mean U_rotor: "+x.get_mean_U_rotor().round(2)+" m/s");
		output.push("OMEGA: "+x.omega.round(1)+" rad/s");
		output.push("LAMBDA: "+x.lambda.round(1));
		output.push("R_blade: "+x.rBlade.round(2)+" m");
		output.push("-----------");
		output.push("THRUST: "+this.thrust.round(2)+" Newtons");
		output.push("POWER : "+this.power.round(2)+" Watts");
		output.push("TORQUE: "+this.torque.round(2)+" Nm");
		//output.push("Cf    : "+this.Cf.round(4));
		if(wakeConvergence.newWakeResults)output.push(wakeConvergence.newWakeResults);
		if(wakeConvergence.isConverged)output.push("WAKE CONVERGED!  STOPPING.");
		output.push("-----------");	
		
		this.resultsOutput(output);
	}
	
	calculateBladePanelSurfaceArea(x){
		//surface area of each quadrilateral blade panel.
		//Panel with nodes ABCD has two triangular subpanels ABC, ACD.
		//Area of triangle with vector sides is:  0.5*|(A-B)x(A-C)|.
		
		var B = x.numBlades * x.nPanels * x.mPanels;

		var node, area;
		for(var i=0; i<B; i++){
			node = x.bladePanel[i].node;
			area = 0.5*node[0].minus(node[1]).cross(node[0].minus(node[2])).NORM();//triangle A.
			area += 0.5*node[0].minus(node[2]).cross(node[0].minus(node[3])).NORM();//triangle B.		
			
			x.bladePanel[i].panel_area = area;
			
			//The original Morino paper multiplied panel covariant vector lengths to get quarter-panel area.  
			//The results are very similar to the above calculation.
			//var a = 4.0*x.bladePanel[i].a1.NORM()*x.bladePanel[i].a2.NORM();		
			//console.log(a+" | "+area+" | "+(a/area));
		}	
	}

	get_Chk_and_Bhk(x, isChk){
		//input: BladeMesh x, boolean isChk
		//output: D_matrix
		
		var t1 = performance.now();
		
		var B = x.nPanels * x.mPanels;

		//The system of blades/wakes is symmetrical.  Because the {phi} vector and the {dphi/dn}
		//vector are repeated, it makes sense to only declare B*B influence matrices, rather
		//than (x.numBlades*B)*(x.numBlades*B) matrices.  This is accomplished by finding the influences over the
		//key blade and then summing them for other blades
		
		//form influence matrices (blade panels)
		var matrix = new D_matrix(B,B);
		var type="Source - ";
		if(isChk)type="Doublet - ";
		
		var k, k_blade, h;
		for(k=0; k<B; k++){
			if(k==0 || (k+1)%100 == 0){
				this.resultsOutput((type+"BLADE k: "+((k+1)*x.numBlades)+" of "+(x.numBlades*B)), true);
			}	

			for(k_blade=0; k_blade<x.numBlades;k_blade++){	
				for(h=0; h<B; h++){
					//doublet influence coefficients for blade - NEWMAN eq 2.14
					if(isChk){
						//matrix is Chk
						if(h==k && k_blade==0)matrix.A[h][k] = 0.5;
						else matrix.A[h][k] += this.doublet(x.bladePanel[k+(k_blade*B)], x.bladePanel[h].P[0]);
					}    	 
					else{
						//matrix is Bhk
						//source influence coefficients for blade - NEWMAN eq 3.10
						matrix.A[h][k] += this.source(x.bladePanel[k+(k_blade*B)], x.bladePanel[h].P[0]);
					}    	
				}	
			}
		}	
		var t2 = performance.now();
		//this.resultsOutput("\n---\ntotal time for Chk and Bhk = "+(t2 - t1)+" milliseconds\n---\n");

		return matrix;
	}
	
	get_Thm(x){
		//input: BladeMesh x
		//output: D_matrix
			
		var B = x.nPanels * x.mPanels;
		var W = x.mPanels * x.turns * x.wPanels;
		
		var Thm = new D_matrix(B,B);
		
		var m, mBlade, h, Fhm, k;
		for(m=0; m<W; m++){
			if(m==0 || (m+1)%100 == 0){
				this.resultsOutput(("Doublet - WAKE  m: "+((m+1)*x.numBlades)+" of "+(W*x.numBlades)), true);
			}
			for(mBlade=0; mBlade<x.numBlades;mBlade++){	
				for(h=0; h<B; h++){
					//doublet influence coefficient for wake - NEWMAN eq 2.14
					Fhm = this.doublet(x.wakePanel[m+(mBlade*W)], x.bladePanel[h].P[0]);

					//and use this to alter the Chk matrix
					//this accounts for the simple steady-state wake influence
					k = x.wakePanel[m].upper;
					Thm.A[h][k] += Fhm;
					k = x.wakePanel[m].lower;
					Thm.A[h][k] -= Fhm;
				}
			}
		}

		return Thm;
	}

	get_local_coordinates(k, Ph){
		//input: SurfacePanel k, Point_3D Ph
		//output Point_3D[]
		
		//a routine for determining the panel-local coordinates of the panel's
		//nodes and the associated field point

		//1: find vector pointing from the centroid of panel k to Ph
		var Pkh = Ph.minus(k.P[0]);

		//2: find unit-direction vectors at panel centroid
		//-->these have already been calculated in blade_panel_outward_normals()
		//   and wake_panel_outward_normals()
		//unit normal is given by k.n

		//3:  find coordinates of Ph in Local coordinate system
		var local_Ph = new Point_3D();
		local_Ph.x = Pkh.dot(k.unit_Xi);
		local_Ph.y = Pkh.dot(k.unit_Eta);
		local_Ph.z = Pkh.dot(k.n);

		//4: find local coordinates of four corner points (local origin at panel centroid)
		var g, local_k = new Array(4); //Point_3D[]
		for(var i=0; i<4; i++){
			g = k.node[i].minus(k.P[0]);
			local_k[i] = new Point_3D();
			local_k[i].x = g.dot(k.unit_Xi);
			local_k[i].y = g.dot(k.unit_Eta);
			local_k[i].z = g.dot(k.n);
		}
		
		//return an array containing the four local corner points, along with the local Ph
		var x = new Array(5); //Point_3D[]
		for(i=0;i<4;i++)x[i] = local_k[i].clone();
		x[4] = local_Ph.clone();

		return x;
	}
	
	doublet(k, Ph){
		//input: SurfacePanel k, Point_3D Ph
		//output: float

		//Doublet strength at a point induced by a quadrilateral panel
		//Newman(1985)
		
		var local_k = new Array(4); //Point_3D[] 
		var localPoints = this.get_local_coordinates(k,Ph); //Point_3D[]
		//unload the results array
		for(var i=0;i<4;i++)local_k[i] = localPoints[i].clone();
		var local_Ph = localPoints[4].clone();
		
		var phi = 0.0, j, delta_Xi, delta_Eta, R, temp;
		for(var i=0; i<4; i++){
			j=0;	
			if(i<3)j=i+1;

			delta_Xi = local_k[j].x - local_k[i].x;
			delta_Eta = local_k[j].y - local_k[i].y;

			R = local_Ph.minus(local_k[i]).NORM();
			temp = delta_Eta*(Math.pow(local_Ph.x - local_k[i].x,2) + local_Ph.z*local_Ph.z);
			temp = temp - delta_Xi*(local_Ph.x - local_k[i].x)*(local_Ph.y - local_k[i].y);
			temp = temp/(R*local_Ph.z*delta_Xi);
			phi = phi + Math.atan(temp);

			R = local_Ph.minus(local_k[j]).NORM();
			temp = delta_Eta*(Math.pow(local_Ph.x - local_k[j].x,2) + local_Ph.z*local_Ph.z);
			temp = temp - delta_Xi*(local_Ph.x - local_k[j].x)*(local_Ph.y - local_k[j].y);
			temp = temp/(R*local_Ph.z*delta_Xi);
			phi = phi - Math.atan(temp);
		}	

		phi = phi/(-4.*Math.PI);

		return phi;
	}
	
	source(k, Ph){
		//input: SurfacePanel k, Point_3D Ph
		//output: float

		//Source strength at a point induced by a quadrilateral panel
		//Newman(1985)
		
		var phi = this.doublet(k,Ph);

		var local_k = new Array(4); //Point_3D[]
		var localPoints = this.get_local_coordinates(k,Ph); //Point_3D[]
		//unload the results array
		for(var i=0;i<4;i++)local_k[i] = localPoints[i].clone();
		var local_Ph = localPoints[4].clone();

		var tol = 1.0e-7;//a small tolerance
		var psi = 0.0, j, A, B, delta_Xi, delta_Eta, R_A, R_B, theta, S_A, temp;
		for(i=0; i<4; i++){
			j=0;
			if(i<3)j=i+1;

			A = local_k[i];
			B = local_k[j];

			delta_Xi = B.x - A.x;
			delta_Eta = B.y - A.y;

			R_A = local_Ph.minus(A).NORM();
			R_B = local_Ph.minus(B).NORM();

			theta=0.;
			if(Math.abs(delta_Xi) < tol){
				if(delta_Eta > 0.0)theta = Math.PI/2.;
				else theta = -Math.PI/2.;
			}
			else {
				theta = Math.atan(delta_Eta/delta_Xi);
				//make sure other quadrants are correct for theta
				 if(delta_Xi < 0.0 && delta_Eta > 0.0)theta = theta + Math.PI;
				 if(delta_Xi < 0.0 && delta_Eta < 0.0)theta = theta - Math.PI;
			}
			S_A = B.minus(A).NORM();

			temp = (local_Ph.x-A.x)*Math.sin(theta) - (local_Ph.y-A.y)*Math.cos(theta);
			psi += temp*Math.log((R_A+R_B+S_A)/(R_A+R_B-S_A));
		}	

		psi = (psi/(-4.*Math.PI)) - local_Ph.z*phi;

		return psi;
	}
	
	boundary_conditions(x){
		//input: BladeMesh x
		//output BladeMesh

		//determine the no-flowthrough velocity boundary condition
		//use eqs 2,3,4,5,6 in Morino et. al. 1980
		
		var B = x.numBlades * x.nPanels * x.mPanels;

		var V_inf = new Point_3D();
		V_inf.x = x.U_0; //the wind blows in the direction of the global X axis
		V_inf.y = 0.0;
		V_inf.z = 0.0; 

		var omega = new Point_3D();
		omega.x = x.omega; //the turbine rotates about the global X axis
		omega.y = 0.0;
		omega.z = 0.0;	

		//find the normal velocity component on each panel collocation point
		for(var i=0; i<B; i++){
			//NOTE:  V_r (in Preuss) will be zero for a rotating frame of reference:
			//       ie:  the blade is stationary and the wind moves past it
			//V_r = CROSS(omega, x.bladePanel[i].P[0]);

			//find velocity of unpeturbed wind in rotating frame of reference
			var V_w = V_inf.minus(omega.cross(x.bladePanel[i].P[0]));
			
			//copy this to the blade panel structure for future use
			x.bladePanel[i].V_w = V_w.clone();

			//determine the component of the inflow normal to this panel
			//the normal potential derivative is set equal (and opposite) to this
			//x.bladePanel[i].dphi_dn = -V_w*x.bladePanel[i].n;
			x.bladePanel[i].dphi_dn = V_w.negative().dot(x.bladePanel[i].n);
		}	
		return x;
	}
	
	get_elemental_quantities(x){
		//input: BladeMesh x
		//output: BladeMesh

		//calculate chord and twist for each blade element
		var k=0, i, blade, TE, LE, b, r, V_w;
		for(blade=0; blade<x.numBlades; blade++){
			for(i=0; i<x.mPanels; i++){
				TE = x.bladePanel[k].node[3];
				if(x.nPanels%2!=0)LE = x.bladePanel[k+(x.nPanels+1)/2].node[3];//odd x.nPanels			
				else LE = x.bladePanel[k+x.nPanels/2].node[3];//even x.nPanels

				//find the vector that points from the TE to the LE
				b = LE.minus(TE);

				//store the chord length and twist of this blade element for future use
				if(blade==0){
					x.element[i].chord = b.NORM();
					r = Math.abs(x.bladePanel[i*x.nPanels].P[0].y);
					V_w = Math.sqrt((r*x.omega)*(r*x.omega) + x.U_rotor[i]*x.U_rotor[i]);
					
					x.element[i].Re = V_w*x.element[i].chord/this.visc;

					if(x.isHelicalWake){
						x.element[i].theta = this.degrees(Math.atan(x.U_rotor[i]/(r*x.omega)));
						x.element[i].twist = this.degrees(Math.atan(-b.z/b.x));
						if(x.element[i].twist < 0.)x.element[i].twist=90.+Math.abs(x.element[i].twist+90.);
						x.element[i].twist = 90.-x.element[i].twist;
						x.element[i].alpha = x.element[i].theta - x.element[i].twist;
					}
					else{
						x.element[i].twist = this.degrees(Math.atan(-b.z/b.x));
						x.element[i].alpha = x.element[i].twist;	    		
						x.element[i].theta=0.0;
					}	    		
				}
				else{
					x.element[blade*x.mPanels + i].chord = x.element[i].chord;
					x.element[blade*x.mPanels + i].twist = x.element[i].twist;
					x.element[blade*x.mPanels + i].alpha = x.element[i].alpha;	    	
					x.element[blade*x.mPanels + i].theta = x.element[i].theta;	    	
					x.element[blade*x.mPanels + i].Re    = x.element[i].Re;	    	
				}

				k += x.nPanels;	    			
			}
		}
		return x;
	}
	
	radians(degs){
		return degs*Math.PI/180.;
	}
	
	degrees(rads){	
		return rads*180./Math.PI;
	}
	
	viscous_boundary_correction(x){
		//input: BladeMesh x
		//output: BladeMesh
		
		//correct the boundary condition in an empirical way for viscous effects

		var RE_CONST1 = 0.2;
		var RE_CONST2 = 1.0;
		var ALPHA_CONST1 = 0.00075;
		var ALPHA_CONST2 = 2.0;
		var F_CONST1 = 1.0;
		var F_CONST2 = 0.0;

	/*	/*------VALUES USED IN OPTIMISATION
			double RE_CONST1=0.3;
			double RE_CONST2=9.0;
			double ALPHA_CONST1=0.0025;
			double ALPHA_CONST2=2.0;
			double F_CONST1=0.2;
			double F_CONST2=1.0;
	*/		
		/*-----ALPHA^2 only correction
			RE_CONST1=0.0;
			RE_CONST2=0.0;
			ALPHA_CONST1=0.95E-3;
			ALPHA_CONST2=2.0;
			F_CONST1=1.0;
			F_CONST2=0.0;
		*/	
		
		//calculate elemental quanities such as chord, twist, angle of attack and Re
		//get_elemental_quantities(X); ->this is now called immediately before current routine
		
		//adjust the boundary condition based on the Reynolds number and the angle of attack
		var k=0, j, RE_MAX=1.E6, blade, element, f_alpha, f_Re, func;
		for(blade=0; blade<x.numBlades; blade++){
			for(element=0; element<x.mPanels; element++){
				f_alpha=ALPHA_CONST1*Math.pow(Math.abs(x.element[element].alpha),ALPHA_CONST2);
			
				f_Re=RE_CONST1*Math.exp(-RE_CONST2*x.element[element].Re/RE_MAX);
				func=F_CONST1*(f_alpha+f_Re) + F_CONST2*(f_alpha*f_Re);		
				if(func<0.0)func=0.0;
				if(func>0.95)func=0.95;
				if(x.element[element].alpha<-5.0)func=0.99;
				for(j=0; j<x.nPanels; j++){
					x.bladePanel[k].dphi_dn *= (1.-func);		
					k++;
				}
			}//end of 'element' loop
		}//end of 'blade' loop	
			
		return x;
	}
	
	get_pressure(x){
		//input: BladeMesh x
		//output: BladeMesh
		
		var B = x.nPanels * x.mPanels;
		var W = x.mPanels * x.turns * x.wPanels;

		//determine the chord-wise location of each blade panel
		x = this.get_x_on_C(x);

		var blade, i, j, k, grad_phi, dist, V_w;
		for(blade=0; blade<x.numBlades; blade++){
			k = blade*B;
			grad_phi = new Point_3D();
			for(i=0; i<x.mPanels; i++){
				for(j=0; j<x.nPanels; j++){
					//find potential derivatives in panel-local coordinates
					//Uses central differences (see Katz and Plotkin eq 12.38, figure 12.25)
					//normal derivative is equal to boundary condition
					grad_phi.z = x.bladePanel[k].dphi_dn;
					//Direction Xi
					dist=0.;
					if(j!=x.nPanels-1 && j!=0){
						dist = x.bladePanel[k+1].a1.NORM() + x.bladePanel[k-1].a1.NORM() + 2.*x.bladePanel[k].a1.NORM();
						grad_phi.x = (x.bladePanel[k+1].phi - x.bladePanel[k-1].phi)/dist;
					}
					else if(j==0){
						dist = x.bladePanel[k].a1.NORM() + x.bladePanel[k+1].a1.NORM();
						grad_phi.x = (x.bladePanel[k+1].phi - x.bladePanel[k].phi)/dist;
					}
					else if(j==x.nPanels-1){
						dist = x.bladePanel[k].a1.NORM() + x.bladePanel[k-1].a1.NORM();
						grad_phi.x = (x.bladePanel[k].phi - x.bladePanel[k-1].phi)/dist;
					}
		
					//Direction Eta
					if(x.mPanels == 1)grad_phi.y = 0.0;
					else{
						if(i!=x.mPanels-1 && i!=0){
							dist = x.bladePanel[k-x.nPanels].a2.NORM() + x.bladePanel[k+x.nPanels].a2.NORM() + 2.*x.bladePanel[k].a2.NORM();
							grad_phi.y = (x.bladePanel[k+x.nPanels].phi - x.bladePanel[k-x.nPanels].phi)/dist;
						}
						else if(i==0){
							dist = x.bladePanel[k+x.nPanels].a2.NORM() + x.bladePanel[k].a2.NORM();
							grad_phi.y = (x.bladePanel[k+x.nPanels].phi - x.bladePanel[k].phi)/dist;
						}
						else if(i==x.mPanels-1){
							dist = x.bladePanel[k].a2.NORM() + x.bladePanel[k-x.nPanels].a2.NORM();
							grad_phi.y = (x.bladePanel[k].phi - x.bladePanel[k-x.nPanels].phi)/dist;
						}
					}

					//resolve the freestream velocity into panel-local components
					V_w = new Point_3D();
					V_w.x = x.bladePanel[k].V_w.dot(x.bladePanel[k].unit_Xi);
					V_w.y = x.bladePanel[k].V_w.dot(x.bladePanel[k].unit_Eta);
					V_w.z = x.bladePanel[k].V_w.dot(x.bladePanel[k].n);

					//determine (steady) pressure coefficient - Katz and Plotkin eq 13.28a, pg430
					x.bladePanel[k].Cp = (1./V_w.dot(V_w)) * (-grad_phi.dot(grad_phi) + 2.*grad_phi.dot(V_w.negative()));
					x.bladePanel[k].pressure = x.bladePanel[k].Cp*0.5*this.rho*V_w.dot(V_w);
		
					k++;
				}
			}
		}

		//make tip and hub pressures equal to the adjacent panels
		for(blade=0; blade<x.numBlades; blade++){
			k = blade*B;
			//hubs
			for(i=0;i<x.nPanels;i++){
				x.bladePanel[i+k].Cp = x.bladePanel[i+k+x.nPanels].Cp;
				x.bladePanel[i+k].pressure = x.bladePanel[i+k+x.nPanels].pressure;	
			}
			//tips
			for(i=B-x.nPanels;i<B;i++){
				x.bladePanel[i+k].Cp = x.bladePanel[i+k-x.nPanels].Cp;
				x.bladePanel[i+k].pressure = x.bladePanel[i+k-x.nPanels].pressure;	
			}
		}

		return x;
	}
	
	get_x_on_C(x){
		//input: BladeMesh x
		//output: BladeMesh
			
		//determine the x/C value of all the panels on the key blade (blade 0)
		
		var k=0, i, j, blade, TE, LE, b, s, P, t, a;
		for(blade=0; blade<x.numBlades; blade++){
			for(i=0; i<x.mPanels; i++){
				TE = x.bladePanel[k].node[3];
				if(x.nPanels%2!=0)LE = x.bladePanel[k+(x.nPanels+1)/2].node[3];//odd x.nPanels			
				else LE = x.bladePanel[k+x.nPanels/2].node[3];//even x.nPanels

				//find the vector that points from the TE to the LE
				b = LE.minus(TE);
				//find the unit vector that points from the TE to the LE
				s = b.divide(b.NORM());
		
				for(j=0; j<x.nPanels; j++){
					//find the midpoint of the current panel's [3]->[0] edge.
					//this approximates the centroid of the panel in 2D
					P = x.bladePanel[k].node[3].plus(x.bladePanel[k].node[0].minus(x.bladePanel[k].node[3]).divide(2));
					//find the vector from the TE to point P
					t = P.minus(TE);
			
					//find the component of t in the direction of the chord
					a = t.dot(s);
			
					//determine x/C from this
					x.bladePanel[k].x_on_C = 1.0 - a/b.NORM();

					k++;
				}
			}
		}
		return x;
	}
	
	get_dP(x){
		//input: BladeMesh x
		//output: D_matrix (column vector)

		var dP = new D_matrix(x.mPanels,1);		
		var upper, lower;
		for(var i=0; i<x.mPanels; i++){
			upper = i*x.nPanels;
			lower = upper + x.nPanels-1;

			dP.A[i][0] = x.bladePanel[upper].Cp - x.bladePanel[lower].Cp;
		}
		
		return dP;
	}
	
	iterate_Thm(x, Thm, dP, iteration){
		//input: BladeMesh x, D_matrix Thm, D_matrix dP, int iteration
		var tol=0.01;
		
		var B = x.nPanels*x.mPanels;	
		
		var j, upper, lower, mult, beta;
		for(var i=0; i<x.mPanels; i++){
			upper = i*x.nPanels;
			lower = upper + x.nPanels-1;

			mult = 1.0 + 0.01/iteration;
			if(Math.abs(dP.A[i][0])>0.5)mult = 1. + 0.03/iteration;	
		
			//determine direction of tweak
			beta=1.0;
			if(dP.A[i][0]<-tol)beta=1./mult;
			else if(dP.A[i][0]>tol)beta=mult;	

			for(j=0; j<B; j++){
				Thm.A[j][upper] *= beta;
				Thm.A[j][lower] *= beta;           	
			}		
		}
				
		return Thm;
	}
	
	calculatePower(x){
		//input: BladeMesh x
			
		//find the power coefficient of this blade via eq 37 in Preuss
		
		var B = x.numBlades * x.nPanels * x.mPanels;

		//define the unit-vector along the X axis (axis of rotation)
		var _i_ = new Point_3D();
		_i_.x = 1.0;
		_i_.y = 0.0;
		_i_.z = 0.0;

		//free-stream dynamic pressure
		var q_dyn = 0.5*this.rho*x.U_0*x.U_0;

		this.torque = 0.0;
		this.thrust = 0.0;
		var force, distance;
		for(var i=0; i<B; i++){
			x.bladePanel[i].force = x.bladePanel[i].pressure * x.bladePanel[i].panel_area;

			force = x.bladePanel[i].n.multiply(-x.bladePanel[i].force);
			distance = _i_.cross(x.bladePanel[i].P[0]);
			
			this.torque += force.dot(distance);
			this.thrust += force.x;
		}
		this.power = this.torque*x.omega;

		var rotor_area = Math.PI*x.rBlade*x.rBlade;

		//power coefficient
		this.Cp = this.power/(q_dyn*rotor_area*x.U_0);

		//thrust coefficient
		this.Cf = this.thrust/(q_dyn*rotor_area);

		return;
	}

	/*
	PotentialFlow.prototype.calculate_delta_phi = function calculate_delta_phi(x){
		//input: BladeMesh x
		//ouput: float
			
		//find the maximum value of delta_phi on the blade for wake iteration
		var delta_phi=0.0, upper, lower, dphi;
		for(var i=0; i<x.mPanels; i++){
		   upper = x.wakePanel[i].upper;
		   lower = x.wakePanel[i].lower;
		   dphi = (x.bladePanel[lower].phi - x.bladePanel[upper].phi)/(x.rBlade*x.U_0);	
		   //use the highest delta_phi value (past 0.5 spanwise) for our U_rotor calculation
		   if(x.bladePanel[upper].r_on_R > 0.2 && dphi>delta_phi)delta_phi = dphi;
		}

		return delta_phi;
	}
	*/
	/*
	PotentialFlow.prototype.calculate_delta_phi = function calculate_delta_phi(x){
		//input: BladeMesh x
		//ouput: float
			
		//find the average value of delta_phi on the blade for wake iteration, weighted for length of TE.
		var upper, lower, dphi;
		var length, total_length = 0., mean_delta_phi = 0.;
		for(var i=0; i<x.mPanels; i++){
		   upper = x.wakePanel[i].upper;
		   lower = x.wakePanel[i].lower;
		   dphi = (x.bladePanel[lower].phi - x.bladePanel[upper].phi)/(x.rBlade*x.U_0);	
		   //average delta_phi, weighted by length of each trailing edge blade panel.
		   if(x.bladePanel[upper].r_on_R > 0.2){
				length = x.bladePanel[lower].node[3].minus(x.bladePanel[lower].node[2]).NORM();	   
				mean_delta_phi += dphi*length;
				total_length += length;
			}	
		}

		return mean_delta_phi/total_length;
	}
	*/
	calculate_delta_phi(x){
		//input: BladeMesh x
		//output: float[]
			
		//find the circulation around each spanwise blade elemement;
		var upper, lower, dphi = new Array(x.mPanels);
		for(var i=0; i<x.mPanels; i++){
		   upper = x.wakePanel[i].upper;
		   lower = x.wakePanel[i].lower;
		   dphi[i] = (x.bladePanel[lower].phi - x.bladePanel[upper].phi)/(x.rBlade*x.U_0);	
		}

		return dphi;
	}

	find_U_rotor(x){
		//input: BladeMesh x
		//output: BladeMesh x	

		//29/3/2001: David's simple wake model
		var delta_phi = this.calculate_delta_phi(x);//find the circulation around the blade.
		var fac, U_inf;
		for(var i=0;i<delta_phi.length;i++){
			fac = x.numBlades*delta_phi[i]*x.lambda/Math.PI;
			U_inf = x.U_0*Math.sqrt(Math.abs(1. - fac));
			if(fac > 1.)U_inf = -U_inf;
			x.U_rotor[i] = (x.U_0 + U_inf)/2.;
		}
				
		return x;		
	}
	
	calculateLiftAndDrag(x){
		//input: BladeMesh x
		//output: BladeMesh
			
		//find the lift and drag (and coefficients) for each blade element
		
		/*
		outfile=fopen("aero_forces.dat","w");
		fresultsOutput.setText(outfile,"R_blade=%2.3lf (metres)\n",R_blade);
		fresultsOutput.setText(outfile,"U_0=%3.3lf (m/s)\n",U_0);
		fresultsOutput.setText(outfile,"LAMBDA=%3.3lf\n",LAMBDA);
		fresultsOutput.setText(outfile,"Blades=%d\n",x.numBlades);	
		fresultsOutput.setText(outfile,"x.nPanels=%d\n",x.nPanels);
		fresultsOutput.setText(outfile,"x.mPanels=%d\n",x.mPanels);	
		fresultsOutput.setText(outfile,"------------------------------------------------------------------\n");
		fresultsOutput.setText(outfile,"M_PANEL\t\ty(m)\t\tf.x(N)\t\tf.y(N)\t\tf.z(N)\n");
		fresultsOutput.setText(outfile,"------------------------------------------------------------------\n");	
		*/
			
		var k=0, j, span, theta, q_dyn, r, V_w, force;
		for(var i=0; i<x.mPanels; i++){
			span = x.bladePanel[i*x.nPanels].node[1].z - x.bladePanel[i*x.nPanels].node[0].z;
			theta=0.;
			q_dyn=0.;
			if(x.isHelicalWake){
				r=x.bladePanel[i*x.nPanels].P[0].z;
				V_w = Math.sqrt((r*x.omega)*(r*x.omega) + x.U_rotor[i]*x.U_rotor[i]);
				theta = this.radians(x.element[i].theta);
				q_dyn = 0.5*this.rho*V_w*V_w*x.element[i].chord*span;
				//this.resultsOutput("["+i+"]:chord="+x.element[i].chord+"  span="+span+"  r="+r+"\n");
			}	  	
			else{
				theta = this.radians(x.element[i].theta);
				q_dyn = 0.5*this.rho*x.U_0*x.U_0*x.element[i].chord*span;
				this.resultsOutput("CALCULATED LIFT AND DRAG FOR STRAIGHT WAKE\n");	  	
			}	  	
		  
			force = new Point_3D();
			for(j=0; j<x.nPanels; j++){
				force = force.minus(x.bladePanel[k].n.multiply(x.bladePanel[k].force));
				k++;
			}
		  
			x.element[i].lift =  force.x*Math.cos(theta) - force.y*Math.sin(theta);
			x.element[i].drag =  force.x*Math.sin(theta) + force.y*Math.cos(theta);

			x.element[i].Cl = x.element[i].lift/q_dyn;
			x.element[i].Cd = x.element[i].drag/q_dyn;
		}

		return x;
	}
	
	wake_convergence(x, U_rotor_old, wake_iteration){
		//input: BladeMesh x, double U_rotor_old, int wake_iteration
		//ouput: object containing results.

		//find the worst ratio between new and old U_rotor iterations.
		var U_rotor_ratio = 1., ratio;
		for(var i=0;i<x.mPanels;i++){
			ratio = x.U_rotor[i]/U_rotor_old.U_minus[i];
			if(x.U_rotor[i] > U_rotor_old.U_minus[i])ratio = U_rotor_old.U_minus[i]/x.U_rotor[i];
			
			if(ratio<U_rotor_ratio)U_rotor_ratio = ratio;
		}	
			
		var oldWakeMessage = "U_rotor convergence = "+(U_rotor_ratio*100.).round(2)+"%";
		var resultsOutput = {isConverged: false, oldWakeResults: oldWakeMessage, newWakeResults: undefined};
		
		if(U_rotor_ratio > 0.999 || (wake_iteration>=5 && U_rotor_ratio > 0.99)){
			resultsOutput.isConverged = true;
			return resultsOutput;//wake has converged
		}	    	
		else{
			//average the new and old U_rotors for this iteration
			for(var i=0;i<x.mPanels;i++){		
				x.U_rotor[i] = (x.U_rotor[i]+U_rotor_old.U_minus[i])/2.;
			}
		}
		//find the ratio
		U_rotor_ratio = 1.;
		for(var i=0;i<x.mPanels;i++){
			ratio = x.U_rotor[i]/U_rotor_old.U_minus[i];
			if(x.U_rotor[i] > U_rotor_old.U_minus[i])ratio = U_rotor_old.U_minus[i]/x.U_rotor[i];
			
			if(ratio<U_rotor_ratio)U_rotor_ratio = ratio;
		}	
		resultsOutput.newWakeResults = "*new* U_rotor convergence = "+(U_rotor_ratio*100.0).round(2)+"%";	
		
		return resultsOutput;//wake has not yet converged
	}
	
	resultsOutput(text, isWipe, isPermanent){
		if(this.isWorker)self.postMessage({'message': text, 'isWipe': isWipe, 'isPermanent': isPermanent});
		else{
			//if(isWipe)document.body.innerHTML = "";
			document.write(text+"<p>");
			//console.log(text);
		}
	}  
	
}//end of PotentialFlow
////
Number.prototype.round = function(places) {
  var result = +(Math.round(this + "e+" + places)  + "e-" + places);
  if(isNaN(result))return 0.;
  return result;
}	
//
//
self.onmessage = function(e){
	if(e.data.bladeMesh)
		var solver = new PotentialFlow(e.data.bladeMesh, true);
}
