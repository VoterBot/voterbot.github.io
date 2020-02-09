class Point_2D{
	constructor(a, b){
		this.x = a;
		this.y = b;
	}

	//methods
	//copy constructor
	clone(a){
		if(a === undefined)a = this;
		return(new Point_2D(a.x, a.y));
	}
	
	//vector subtraction
	minus(p){
		return new Point_2D(this.x - p.x, this.y - p.y);
	}
	
	//vector addition
	plus(p){
		return new Point_2D(this.x + p.x, this.y + p.y);
	}
	
	//multiply by a constant
	multiply(a){
		return new Point_2D(this.x*a, this.y*a);
	}
	
	//divide by a constant
	divide(a){
		return new Point_2D(this.x/a, this.y/a);
	}
	
	//dot product
	dot(p){
		return (this.x*p.x + this.y*p.y);
	}
	
	//Norm (length of vector)
	NORM(){
		return Math.sqrt(this.x*this.x + this.y*this.y);
	}
	
	//distance between two points
	distance2D(p){
		return(Math.sqrt(Math.pow(this.x - p.x,2) + Math.pow(this.y - p.y,2)));
	}
	
	//a string with a label for output
	toString(label){
		return (label+"=("+this.x+","+this.y+")");
	}
	output(label){
		alert(this.toString(label));	
	}
}	
