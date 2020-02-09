"use strict";
class LeastSquaresCurve{
	constructor(Q, n, p){
		//least-squares curve fitting to point data. A9.6 p417 P&T.
		//input: Q, the point data to fit a curve to: Q[i]: 2D/3D point and weight, stored as WeightedPoint_2D/3D; 
		//											  Q[i].d: optional derivative at point & weight.
		//input: n, the high index of control point vector for the new NURBS curve; p, the degree of the curve.
		//Output:  A **non-rational** NURBS curve, which means a Bspline curve.  

		//initialise the routine
		if(Q)this.c = this.createLeastSquaresCurve(Q,n,p); //a new NURBSCurve
	}	
	
	//methods
	createLeastSquaresCurve(Q, n, p){
		var i, j;

		var ru = -1; 
		var rc = -1;
		var su = -1;
		var sc = -1;
		
		for(i=0; i<Q.length; i++){
			if(Q[i].w > 0.0)ru++;
			else rc++;
			
			if(Q[i].d){
				if(Q[i].d.w > 0.0)su++;
				else sc++;
			}	
		}

		var mu = ru+su+1;
		var mc = rc+sc+1;
		
		//is the dataset 2D or 3D?
		var is3D = true;
		if(Q[0].p.z === undefined)is3D = false;
		
		//set up some matrices.  D_matrix is float[][].  P_matrix is WeightedPoint_2D (or _3D)[][]
		var N = new D_matrix(mu+1,n+1);
		var M = new D_matrix(mc+1,n+1);
		var S = new D_matrix(mu+1,1); S.P_matrix(is3D);
		var TT = new D_matrix(mc+1,1); TT.P_matrix(is3D);
		var A = new D_matrix(mc+1,1); A.P_matrix(is3D);
		var W = new D_matrix(mu+1,mu+1);
		var P = new D_matrix(n+1,1); P.P_matrix(is3D);
	  
		var NtWN = new D_matrix(n+1,n+1);
		var NtWS = new D_matrix(n+1,1); NtWS.P_matrix(is3D);
		
		//compute and load parameters uk into ub[] (EQ. 9.5)
		var ub = this.calculateParameters(Q); //float[]
		
		//compute and load the knots into U[] (EQ 9.68 and 9.69)
		var U = this.calculateKnots(n,p,ub); //float[]

		//counters up to mu and mc
		var mu2=0; //int
		var mc2=0; //int
		
		var tempCurve = new NURBSCurve();
		tempCurve.n = n;
		tempCurve.p = p;
		tempCurve.U = U;
		
		var Nip, span, Nd;
		for(i=0; i<Q.length; i++){
			span = tempCurve.findSpan(ub[i]);
			Nd = tempCurve.dersBasisFuns(span, ub[i], 1);
			
			if(Q[i].w > 0.){ //unconstrained point
				W.A[mu2][mu2] = Q[i].w;
				
				//load mu2th row of N[][] with this span's basis functions.			
				for(j=0;j<=p;j++)N.A[mu2][j+span-p] = Nd[0][j]; 

				S.A[mu2][0] = Q[i].pw(); //matrix routines require a weighted point, so use .pw for convenience.
				mu2++;
			}
			else{ //constrained point
				//load mc2th row of M[][] with this span's basis functions.		
				for(j=0;j<=p;j++)M.A[mc2][j+span-p] = Nd[0][j];
				
				TT.A[mc2][0] = Q[i].clone();
				mc2++;
			}

			if(Q[i].d){ //derivative specified.
				if(Q[i].d.w > 0.){ //unconstrained derivative.
					W.A[mu2][mu2] = Q[i].d.w;
					
					//load mu2th row of N[][] with this basis function derivatives.
					for(j=0;j<=p;j++)N.A[mu2][j+span-p] = Nd[1][j];
					
					S.A[mu2][0] = Q[i].d.pw();
					mu2++;
				}
				else{ //constrained derivative.
					//load mc2th row of M[][] with this span's basis function derivatives.			
					for(j=0;j<=p;j++){		
						M.A[mc2][j+span-p] = Nd[1][j]; 
					}
					
					TT.A[mc2][0] = Q[i].d.clone();
					mc2++;
				}
			}
		} // end of for-loop i=0,..,r
		
		//compute NtWN & NtWS matrices
		var Nt = N.T();
		NtWN = (Nt.cross(W)).cross(N);
		NtWS = (Nt.cross(W)).cross(S);
		
		if(mc < 0){ //no constraints
			P = NtWN.I().cross(NtWS);

			var c = new NURBSCurve();
			c.initialise(n,p);
			for(i=0; i<=n; i++){
				c.P[i] = P.A[i][0].clone();
				c.P[i].w = 1.0; //This least-squares routine returns a non-rational curve.
			}	
			c.U = new Array(U.length);
			for(i=0;i<U.length;i++)c.U[i]=U[i];
		  
			return c;
		}
		
		//compute the inverse of NtWN and store this back in NtWN
		NtWN = NtWN.I();
		
		var temp_d = (M.cross(NtWN)).cross(M.T()); //D_matrix
		
		var temp_p = ((M.cross(NtWN)).cross(NtWS)).minus(TT);
		
		A = temp_d.I().cross(temp_p); //eq 9.75
		
		P = NtWN.cross(NtWS.minus(M.T().cross(A))); //eq 9.74
		
		var c = new NURBSCurve();
		c.initialise(n,p);
		for(i=0; i<=n; i++){
			c.P[i] = P.A[i][0].clone();
			c.P[i].w = 1.0; //This least-squares routine returns a non-rational curve.			
		}	
		c.U = new Array(U.length);
		for(i=0;i<U.length;i++)c.U[i]=U[i];
		
		return c;
	}  
	
	calculateParameters(Q){
		//input: Point_2D[]; output: float[]
		//chord length method.  Better than equally-spaced.  Not so good at sharp turns.
		//find parameters based on the chord length (EQ 9.4 and 9.5)
		var d=0.;
		for(var i=1;i<Q.length;i++){
			d+=(Q[i].p.minus(Q[i-1].p)).NORM();
		}
		
		var u = new Array(Q.length);
		u[0]=0.;
		u[u.length-1]=1.;
		  
		for(var i=1;i<u.length-1;i++){
			u[i] = u[i-1] + (Q[i].p.minus(Q[i-1].p)).NORM()/d;
		}  
		  
		return u;
	}

	/*
	calculateParameters(Q){
		//input: Point_2D[]; output: float[]
		//Centripital method.  Better at sharp turns.  Supposedly.
		//find parameters based on the chord length (EQ 9.6)
		var d=0.;
		for(var i=1;i<Q.length;i++){
			d+=Math.sqrt(Q[i].p.minus(Q[i-1].p).NORM());
		}
		
		var u  = new Array(Q.length);
		u[0]=0.;
		u[u.length-1]=1.;
		  
		for(var i=1;i<u.length-1;i++){
			u[i] = u[i-1] + Math.sqrt(Q[i].p.minus(Q[i-1].p).NORM())/d;
		}  
		  
		return u;
	}
	*/
	
	calculateKnots(n, p, ub){
		//input: int n, int p, float[] ub
		//output: float[]
		//EQ 9.68 and 9.69
		
		var m = ub.length-1;//high-index of parameter vector (also high-index of points to be interpolated)
		  
		var d = (m+1.)/(n-p+1.);
		  
		var U = new Array(n+p+2); //float[]
		var i, j, a;
		for(j=1;j<=n-p;j++){
			i = Math.floor(j*d);
			a = j*d - i;
			U[p+j] = (1.-a)*ub[i-1]+a*ub[i];  
		}

		//p+1 repeated knots at the endpoints
		for(i=0;i<=p;i++){
			U[i]=0.;
			U[i+n+1]=1.;
		}  

		return U;
	}
	
	extract2DCurve(type){
		if(this.c.P[0].p.z === undefined)return this.c; //it's already a 2D curve.

		var curve = this.c.clone();
		for(var i=0;i<curve.P.length;i++){
			if(type.isXY)curve.P[i] = new WeightedPoint_2D(this.c.P[i].p.x, this.c.P[i].p.y, this.c.P[i].w);
			else if(type.isXZ)curve.P[i] = new WeightedPoint_2D(this.c.P[i].p.x, this.c.P[i].p.z, this.c.P[i].w);		
		}
		
		return curve;
	}

	clone(a){
		if(a === undefined)a = this;
		var x = new LeastSquaresCurve();
		x.c = NURBSCurve.prototype.clone(a.c);
		return x;
	}	
}