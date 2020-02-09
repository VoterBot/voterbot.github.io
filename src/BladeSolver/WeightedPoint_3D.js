class WeightedPoint_3D{
	constructor(x,y,z,weight){
		this.p = new Point_3D(x,y,z);
		this.w = weight; //float
	}
	
	//methods
	//copy constructor
	clone(a){
		if(a === undefined)a = this;
		return new WeightedPoint_3D(a.p.x, a.p.y, a.p.z, a.w);
	}
	
	//create a weighted point from a normal 3D point
	copy_Point_3D(point){
		return new WeightedPoint_3D(point.x,point.y,point.z,0.);
	}

	//vector addition
	plus(wpoint){
		var result = new WeightedPoint_3D(0.,0.,0.,0.);
		result.p = this.p.plus(wpoint.p);
		result.w = this.w+wpoint.w;
	  
		return result;
	}

	//vector subtraction
	minus(wpoint){
		var result = new WeightedPoint_3D(0.,0.,0.,0.);
		result.p = this.p.minus(wpoint.p);
		result.w = this.w-wpoint.w;
	  
		return result;
	}

	//multiply point and weight by a constant
	multiply(a){
		var result = new WeightedPoint_3D(0.,0.,0.,0.);
		result.p = this.p.multiply(a);
		result.w = this.w*a;
	  
		return result;
	}

	//divide point and weight by a constant
	divide(a){
		var result = new WeightedPoint_3D(0.,0.,0.,0.);
		result.p = this.p.divide(a);
		result.w = this.w/a;
	  
		return result;
	}

	//a function to return a weighted point (the point multiplied by the weight)
	pw(){
		var result = new WeightedPoint_3D(0.,0.,0.,0.);
		result.p = this.p.multiply(this.w);
		result.w = this.w;

		return result;
	}

	//Norm (length of vector)
	NORM(){
		return this.p.NORM();
	}

	toString(name){
		return (name+"=("+this.p.x+","+this.p.y+","+this.p.z+","+this.w+")<br>");
	}

	output(name){
		document.write(this.toString(name));
	}
}   
