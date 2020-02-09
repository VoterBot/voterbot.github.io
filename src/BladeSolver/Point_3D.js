class Point_3D{
	constructor(a, b, c){
		this.x = a;
		this.y = b;
		this.z = c;	
	}	
	
	//methods
	//copy constructor
	clone(a){
		if(a === undefined)a = this; //default copy constructor.
		return(new Point_3D(a.x, a.y, a.z));
	}

	//vector subtraction
	minus(p){
		return new Point_3D(this.x - p.x, this.y - p.y, this.z - p.z);
	}

	//vector addition
	plus(p){
		return new Point_3D(this.x + p.x, this.y + p.y, this.z + p.z);
	}

	//multiply by a constant
	multiply(a){
		return new Point_3D(this.x*a, this.y*a, this.z*a);
	}

	//negative
	negative(){
		return this.clone().multiply(-1.);
	}


	//divide by a constant
	divide(a){
		return new Point_3D(this.x/a, this.y/a, this.z/a);
	}

	//dot product
	dot(p){
		return this.x*p.x + this.y*p.y + this.z*p.z;
	}

	//cross product
	cross(p){
		//see page 351 in "Advanced Engineering Mathematics - Zill and Cullen
		var result = this.clone();

		result.x = this.y*p.z - this.z*p.y;
		result.y = this.z*p.x - this.x*p.z;
		result.z = this.x*p.y - this.y*p.x;

		return result;
	}

	//Norm (length of vector)
	NORM(){
		return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
	}

	//distance between two points
	distance3D(p){
		//return(Math.sqrt(Math.pow(this.x - p.x,2) + Math.pow(this.y - p.y,2)));
		var dx = this.x - p.x;
		var dy = this.y - p.y;
		var dz = this.z - p.z;	
		return Math.sqrt(dx*dx + dy*dy + dz*dz);		
	}

	//a string with a label for output
	toString(label){
		return (label+"=("+this.x+","+this.y+","+this.z+")");
	}
};
