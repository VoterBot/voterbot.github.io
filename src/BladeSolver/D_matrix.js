class D_matrix{
	constructor(rows_, cols_){
		//matrix storage for floating point matrix algebra
		this.rows = rows_;
		this.cols = cols_;	
		this.A = create2DArray(this.rows, this.cols);
	}	
	
	//methods
	//copy constructor
	clone(){
		var x = new D_matrix(this.rows, this.cols);
		
		var i, j;
		for(i=0; i<this.rows; i++)
			for(j=0; j<this.cols; j++)
				x.A[i][j] = this.A[i][j];
				
		x.rows = this.rows;
		x.cols = this.cols;
				
		return x;		
	}
	
	//matrix addition
	plus(x){
		var w = x.clone();

		var i, j;
		for(i=0; i<w.rows; i++){
			for(j=0; j<w.cols; j++){
				w.A[i][j] = w.A[i][j].plus(this.A[i][j]);
			}
		}
		return w;
	}
	
	//matrix subtraction
	minus(x){
		var w = x.clone();

		var i, j;
		for(i=0; i<w.rows; i++){
			for(j=0; j<w.cols; j++){
				w.A[i][j] = this.A[i][j].minus(w.A[i][j]);
			}
		}
		return w;
	}
	
	//scalar multiplication
	multiply(x){
		var w = this.clone();

		var i, j;
		for(i=0; i<w.rows; i++){
			for(j=0; j<w.cols; j++){
				w.A[i][j] = w.A[i][j].multiply(x);
			}
		}
		return w;
	}
	
	//get -X
	negative(){
		var zero = new D_matrix(this.rows, this.cols);
		return zero.minus(this);
	}	
	
	//cross product: this * X
	cross(x){
		var w = new D_matrix(this.rows, x.cols);
		if(x.isP_matrix())w.P_matrix(x.A[0][0].p.z !== undefined); //if we're multiplying a float[][] by a WeightedPoint[][]
		
		var i, j, k;
		for(k=0; k<this.cols; k++){
			for(i=0; i<w.rows; i++){
				for(j=0; j<w.cols; j++){
					//w.A[i][j] += this.A[i][k]*x.A[k][j];
					w.A[i][j] = w.A[i][j].plus(x.A[k][j].multiply(this.A[i][k]));
				}
			}
		}	
		return w;
	}
	
	//find the transpose of this matrix
	T(){
		var w = new D_matrix(this.cols, this.rows);

		var i, j;
		for(i=0; i<this.rows; i++){
			for(j=0; j<this.cols; j++){
				w.A[j][i] = this.A[i][j];
			}
		}		

		return w;
	}
	
	//find the inverse of this matrix
	I(){	
		var a = this.clone();
		
		//decompose the matrix once
		var indx = a.ludcmp(); //int[]

		var y = new D_matrix(this.rows, this.cols);
		//find inverse by columns
		var col  = new Array(this.rows); //float[]
		var i, j;
		for(j=0; j<this.rows; j++){
			for(i=0; i<this.rows; i++)col[i]=0.;
			col[j]=1.;
			col = a.lubksb(indx, col);
			for(i=0; i<this.rows; i++)y.A[i][j]=col[i];
		}	

		return y;	
	}
	
	//LU Decomposition.
	ludcmp(){
		//from Numerical Recipes in C
		/*
		Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition
		of a rowwise permutation of itself. a and n are input. a is output, arranged
		as in equation (2.3.14) above;  indx[0..n-1] is an output vector that records the
		row permutation eected by the partial pivoting; d is output as  1 depending on
		whether the number of row interchanges was even or odd, respectively. This routine
		is used in combination with lubksb to solve linear equations or invert a matrix.
		*/
		
		var tiny = 1.E-10;

		var vv = new Array(this.rows); //float[]

		var i, j, big, temp;
		for(i=0; i<this.rows; i++){//loop over rows to get the implicit scaling information
			big=0.;
			for(j=0; j<this.rows; j++){
				temp = Math.abs(this.A[i][j]);
				if (temp > big)big=temp;
			}	
			if(big == 0.){
				console.log("OOPS!  Singular Matrix.  Stopping\n");
				return [];
			}
			vv[i] = 1./big;	//save the scaling
		}
		
		var imax=0; //int
		var dum=0.; //float
		var indx = new Array(this.rows); //int[]
		var k, sum;
		for(j=0; j<this.rows; j++){	//this is the loop over columns of Crout's method.
			for(i=0; i<j; i++){ //eq. 2.3.12 except for 1=j
				sum = this.A[i][j];
				for(k=0; k<i; k++)sum -= this.A[i][k]*this.A[k][j];
				this.A[i][j]=sum;
			}
			big=0.;	//initialise for the search for largest pivot element
			for(i=j; i<this.rows; i++){          //this is i=j of eq. 2.3.12 and
				sum = this.A[i][j];     //i=j+1...N of equation 2.3.13
				for(k=0; k<j; k++)sum -= this.A[i][k]*this.A[k][j];
				this.A[i][j]=sum;
				if((dum=vv[i]*Math.abs(sum)) >= big){
					//is the figure of merit for the pivot better than
					//the best so far?
					big = dum;
					imax = i;
				}
			}
			if(j != imax){		     //Do we need to interchange rows?
				for(k=0; k<this.rows; k++){ //Yes, so do so...
					dum = this.A[imax][k];
					this.A[imax][k] = this.A[j][k];
					this.A[j][k] = dum;
				}
				vv[imax]=vv[j];	    //Also interchange the scale factor
			}
			indx[j]=imax;
			if(this.A[j][j] == 0.)this.A[j][j] = tiny;
			//if the pivot element is zero the matrix is singular (at least to
			//the precision of the algorithm).  For some applications on singular
			//matrices it is desirable to substitute TINY for zero.
			if(j != this.rows){        //now, finally, divide by the pivot element
				dum = 1./(this.A[j][j]);
				for(i=j+1; i<this.rows; i++)this.A[i][j] *= dum;
			}
		} //go back for the next column in the reduction.

		return indx;
	}
	
	//Solve linear equations by back substitution
	lubksb(indx, b){
		//input: indx (int[]); b (float[])
		//returns: float[]
		//from Numerical Recipes in C
		/*
		Solves the set of n linear equations A  X = B. Here a[0..n-1][0..n-1] is input,
		not as the matrix A but rather as its LU decomposition, determined by the
		routine ludcmp. indx[0..n-1] is input as the permutation vector returned by
		ludcmp. b[0..n-1] is input as the right-hand side vector, and returns with
		the solution vector X. a, n, and indx are not modified by this routine
		and can be left in place for successive calls with different right-hand sides b.
		This routine takes into account the possibility that b will begin with many
		zero elements, so it is efficient for use in matrix inversion.
		*/

		var sum=0.; //float
		var i, j, ii=0, ip; //int
		for(i=1; i<=this.rows; i++){	//when ii is set to a positive value, it will become
			ip = indx[i-1];	//the index of the first nonvanishing element of b.
			sum = b[ip];	//We now do the forward substitution, eq 2.3.6.  The
			b[ip] = b[i-1];	//only new wrinkle is to unscramble the permutation
			if(ii>0)	//as we go
				for(j=ii; j<= i-1; j++)sum -= this.A[i-1][j-1]*b[j-1];
			else if(sum>0.) ii = i; //a nonzero element was encountered, so from
			b[i-1] = sum;		//now on we will have do do the sums in the loop above
		}

		for(i=this.rows-1; i>=0; i--){		//now we do the back-substitution, eq 2.3.7
			sum = b[i];
			for(j=i+1; j<this.rows; j++)sum -= this.A[i][j]*b[j];
			b[i]=sum/this.A[i][i];	//store a component of the solution vector X
		}//all done!

		return b;
	}
	
	//solve a linear set of equations Ax=b
	solve(Y){
		//returns the vector x

		var A = this.clone(); //copy of X
		var b = Y.clone(); //copy of Y		

		var N = A.rows;
		var col = new Array(N); //float[]

		//load the "col" vector with the values in b
		var i;
		for(i=0; i<N; i++)col[i] = b.A[i][0];

		//decompose the matrix once
		var indx = A.ludcmp(); //int[]
		//and back substitute our RHS vector
		col = A.lubksb(indx, col);
		
		//copy the result to our solution vector
		var x = new D_matrix(N,1);
		for(i=0; i<N; i++)x.A[i][0] = col[i];
		
		return x;
	}	
	
	//convert current float[][] matrix into a Point_2D (or _3D) [][] matrix.	
	P_matrix(is3D){
		var p;
		if(is3D)p = new WeightedPoint_3D(0.,0.,0.,1.);
		else p = new WeightedPoint_2D(0.,0.,1.);	
		
		var i, j;
		for(i=0; i<this.rows; i++)
			for(j=0; j<this.cols; j++)
				this.A[i][j] = p.clone();
	}
	
	//if the first element is a number, it's a standard floating point/integer matrix.
	isP_matrix(){
		return ((typeof this.A[0][0])!=="number");
	}

	output(title){
		var result = ("<br>"+title+"<br>");
		var i, j, r;
		for(i=0;i<this.rows;i++){
			r = "";
			for(j=0;j<this.cols;j++){
				if(this.isP_matrix())this.A[i][j].output("["+i+"]["+j+"]");
				else r+=(" "+this.A[i][j].round(3));
			}
			if(this.isP_matrix())document.write("<br>");
			result += (r+"<br>");
		}
		console.log(result);
	}
}
//
function create2DArray(a, b){
	var x = new Array(a);
	var i, j;
	for(i=0;i<a;i++){
		x[i] = new Array(b);
		for(j=0;j<b;j++)x[i][j]=0.;
	}	
	
	return x
};
//
Number.prototype.round = function(places) {
  var result = +(Math.round(this + "e+" + places)  + "e-" + places);
  if(isNaN(result))return 0.;
  return result;
};		
//create operators for plus, minus and multiply, and add them to the Number class.  These mirror
//the functions in the Point_2D class, and allows the same routines to be used on floats and Point_2D.
Number.prototype.plus = function(x) {
	return this+x;
};
Number.prototype.minus = function(x) {
	return this-x;
};
Number.prototype.multiply = function(x) {
	return this*x;
};

/*	
function test1(){	
		var A = new D_matrix(3,3);
		
		A.A[0][0]=2;
		A.A[0][1]=0;
		A.A[0][2]=1;
		A.A[1][0]=-2;
		A.A[1][1]=3;
		A.A[1][2]=4;
		A.A[2][0]=-5;
		A.A[2][1]=5;
		A.A[2][2]=6;
		
		A.output("A");	
		
		var B = A.I();
		
		A.output("A");
		B.output("B = I(A)");
		
		var C = A.cross(B);
		
		C.output("C = AxB");
}

function test3(){	
		var A = new D_matrix(3,3);
		
		A.A[0][0]=1;
		A.A[0][1]=2;
		A.A[0][2]=3;
		A.A[1][0]=0;
		A.A[1][1]=4;
		A.A[1][2]=5;
		A.A[2][0]=1;
		A.A[2][1]=0;
		A.A[2][2]=6;
		
		A.output("A");	
		
		var B = A.I().multiply(11);

		
		A.output("A");
		B.output("B = I(A)*11");
		
		var C = A.cross(B);
		
		C.output("C = AxB");
}

function test2(){
		var X = new D_matrix(4,3);
		
		X.A[0][0] = 14;
		X.A[0][1] = 9;
		X.A[0][2] = 3;
		X.A[1][0] = 2;
		X.A[1][1] = 11;
		X.A[1][2] = 15;
		X.A[2][0] = 0;
		X.A[2][1] = 12;
		X.A[2][2] = 17;
		X.A[3][0] = 5;
		X.A[3][1] = 2;
		X.A[3][2] = 3;
		
		var B = new D_matrix(3,2);
		B.A[0][0] = 12;
		B.A[0][1] = 25;
		B.A[1][0] = 9;
		B.A[1][1] = 10;						
		B.A[2][0] = 8;
		B.A[2][1] = 5;
		
		var C = X.cross(B);

		X.output("X");
		B.output("B");
		C.output("C");		
}	

function test4(){
	//solve an example linear system.  From http://en.wikipedia.org/wiki/System_of_linear_equations
	//Ax=b.  Solve for x.
	//A = 	[1 3 -2]
	//		[3 5 6 ]
	//		[2 4 3 ]
	//
	//b = 	[5]
	//		[7]
	//		[8]
	//
	//x = 	[-15]
	//		[8]
	//		[2]
	
	var startTime = window.performance.now();
	
	var A = new D_matrix(3,3);
	
	A.A[0][0] = 1.;
	A.A[0][1] = 3.;
	A.A[0][2] = -2.;		
	A.A[1][0] = 3.;
	A.A[1][1] = 5.;
	A.A[1][2] = 6.;		
	A.A[2][0] = 2.;
	A.A[2][1] = 4.;
	A.A[2][2] = 3.;		
	
	var b = new D_matrix(3,1);
	
	b.A[0][0] = 5.;
	b.A[1][0] = 7.;
	b.A[2][0] = 8.;
	
	A.output("A");
	
	b.output("b");
	
	var x = A.solve(b);
	
	x.output("x");
	
	A.negative().output("-A");
	
	A.I().cross(A).output("A<sup>-1</sup>xA");	
	
	A.T().output("T(A)");
	
	var timeNow = window.performance.now();

	document.write("<br>----<br>calculation took "+(timeNow-startTime)+" milliseconds.");
}	
*/