class WeightedPoint_2D{
	constructor(x,y,weight){
		this.p = new Point_2D(x,y);
		this.w = weight; //float
	}	
	
	//methods
	//copy constructor
	clone(a){
		if(a === undefined)a = this;
		return new WeightedPoint_2D(a.p.x, a.p.y, a.w);
	}

	//create a weighted point from a normal 2D point
	copy_Point_2D(point){
		return new WeightedPoint_2D(point.x,point.y,0.);
	}

	//vector addition
	plus(wpoint){
		var result = new WeightedPoint_2D(0.,0.,0.);
		result.p = this.p.plus(wpoint.p);
		result.w = this.w+wpoint.w;

		return result;
	}  

	//vector subtraction
	minus(wpoint){
		var result = new WeightedPoint_2D(0.,0.,0.);
		result.p = this.p.minus(wpoint.p);
		result.w = this.w-wpoint.w;

		return result;
	}  

	//multiply point and weight by a constant
	multiply(a){
		var result = new WeightedPoint_2D(0.,0.,0.);
		result.p = this.p.multiply(a);
		result.w = this.w*a;

		return result;
	}  

	//divide point and weight by a constant
	divide(a){
		var result = new WeightedPoint_2D(0.,0.,0.);
		result.p = this.p.divide(a);
		result.w = this.w/a;

		return result;
	}

	//distance between two weighted points.  Treat the weight just like a coordinate.
	distance3D(p){
		var xx = this.p.x - p.p.x;
		var yy = this.p.y - p.p.y;
		var ww = this.w - p.w;
		return(Math.sqrt(xx*xx + yy*yy + ww*ww));
	}

	//a function to return a weighted point (the point multiplied by the weight)
	pw(){
		var result = new WeightedPoint_2D(0.,0.,0.);
		result.p = this.p.multiply(this.w);
		result.w = this.w;

		return result;
	}
	
	output(name){
		alert(name+"=("+this.p.x+","+this.p.y+","+this.w+")");
	}  
}
