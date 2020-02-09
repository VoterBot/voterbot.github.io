"use strict";
class PlotPanel{
	constructor(parameters){
		/*
			parameters = {
				plotContainer: div,
				plotWidth: ,
				plotMargin: {left:, right:, top:, bottom:},
				numPlotsAcrossPage:  ,
				hasBorder: boolean,
			}
		*/
		this.parameters = parameters;
		this.currentIndex = 0;

		//container border.
		if(this.parameters.hasBorder){
			this.parameters.plotContainer.style.borderColor = "black";
			this.parameters.plotContainer.style.borderStyle = "solid";
			this.parameters.plotContainer.style.borderWidth = "5px";
		}
		
		//set aspect ratio;
		if(this.parameters.aspectRatio === undefined)this.parameters.aspectRatio = 1.;//square plot.
		this.parameters.plotHeight = this.parameters.plotWidth*this.parameters.aspectRatio;

		//a container for each set of plot data.
		this.plots = {};
	}

	//methods
	//setter.
	addPlot(plotData){
		var plot = new Plot(plotData, this);
		
		//store the data.
		this.plots[plot.name] = plot;
		//draw the plot.
		plot.drawPlot();
	}
}	
////
////
class Plot{
	constructor(data, plotPanel){
		var self = this;
		var p = plotPanel.parameters;

		this.data = data;
		this.plotPanel = plotPanel;
		
		this.name = 'canvas'+plotPanel.currentIndex++;
		
		//create a new plot canvas.
		var canvas = document.createElement('canvas');
		canvas.setAttribute('id', this.name); //give the canvas a unique name.
		//if the plot has a mouse-handling callback defined, attach some mouse event listeners.
		if(data.callback){
			canvas.addEventListener("mousedown", function(e){self.mouseDownListener(e);}, false);
			canvas.addEventListener("mouseup", function(e){self.mouseUpListener(e);}, false);		
			canvas.addEventListener("mousemove", function(e){self.mouseMoveListener(e);}, false);
		}	
		
		//setup.
		this.ctx = canvas.getContext("2d");
		canvas.width = p.plotWidth;
		canvas.height = p.plotHeight;
		
		//add the canvas to the page.
		p.plotContainer.appendChild(canvas);
		
		//resize the container div to accommodate this new plot.
		var numPlotsDownPage = Math.ceil(plotPanel.currentIndex/p.numPlotsAcrossPage);
		p.plotContainer.style.width = p.plotWidth*p.numPlotsAcrossPage+"px";
		p.plotContainer.style.height = p.plotHeight*numPlotsDownPage+"px";
		
		//set up axes and screen coordinates for this data.
		var axes = this.getAxes(data);
		this.drawBoundary = {
			topLeftX: p.plotMargin.left,
			topLeftY: p.plotMargin.top,
			width: p.plotWidth - p.plotMargin.left - p.plotMargin.right,
			height: p.plotHeight - p.plotMargin.top - p.plotMargin.bottom,
			xAxis: axes.xAxis,
			yAxis: axes.yAxis,
		}

		//storage for mouse point.
		this.mousePoint = new Point_2D();
	}

	//methods
	//
	toScreenCoord(coord){
		var d = this.drawBoundary;
		var screenCoord = {};
		if(coord.x != undefined)
			screenCoord.x = d.topLeftX + d.width*(coord.x - d.xAxis.min)/(d.xAxis.max - d.xAxis.min);
		//y-axis needs to be inverted.
		if(coord.y != undefined)
			screenCoord.y = d.topLeftY + d.height - d.height*(coord.y - d.yAxis.min)/(d.yAxis.max - d.yAxis.min);
		return screenCoord;
	}

	fromScreenCoord(screenCoord){
		var d = this.drawBoundary;
		var coord = {};
		if(screenCoord.x != undefined)
			coord.x = d.xAxis.min + ((screenCoord.x - d.topLeftX)/d.width)*(d.xAxis.max - d.xAxis.min);
		//y-axis needs to be inverted.
		if(screenCoord.y != undefined)
			coord.y = d.yAxis.min + ((screenCoord.y - d.topLeftY - d.height)/(-d.height))*(d.yAxis.max - d.yAxis.min);
		return coord;
	}

	getMouseHover(){
		//store the screen locations of control points for bspline curves, if we haven't already.
		if(this.screenControlPoints === undefined && this.data.bsplineCurveFit){
			this.screenControlPoints = [];
			var p;
			for(var i=0;i<this.data.bsplineCurveFit.P.length;i++){
				p = this.toScreenCoord(this.data.bsplineCurveFit.P[i].p);
				this.screenControlPoints.push(new Point_2D(p.x, p.y));
			}
		}
		
		//check to see if the mouse point is close to a control point.
		if(this.screenControlPoints){
			this.hoverIndex = undefined;
			var dist, tol = 5;
			for(var i=0;i<this.screenControlPoints.length;i++){
				dist = this.mousePoint.minus(this.screenControlPoints[i]).NORM();
				if(dist<=tol){
					this.hoverIndex = i;
					break;
				}
			}
		}
	}

	dragControlPoint(allowDrag){
		//console.log("DRAGGING: "+this.draggingIndex+" to:"+this.mousePoint.toSource());
		if(allowDrag === undefined)allowDrag = {x:true, y:true};
		
		if(allowDrag.x)this.screenControlPoints[this.draggingIndex].x = this.mousePoint.x;
		if(allowDrag.y)this.screenControlPoints[this.draggingIndex].y = this.mousePoint.y;
		
		var p = this.fromScreenCoord(this.mousePoint);
		if(allowDrag.x)this.data.bsplineCurveFit.P[this.draggingIndex].p.x = p.x;
		if(allowDrag.y)this.data.bsplineCurveFit.P[this.draggingIndex].p.y = p.y;
	};

	getCanvas(event){
		//a helper for listeners:  extracts canvas, context and plot from mouse event.
		var c = {};
		c.canvas = event.target;
		c.ctx = c.canvas.getContext('2d');
		c.plot = this.plotPanel.plots[c.ctx.canvas.id];
		
		return c;
	}

	mouseDownListener(event){
		//console.log(event);

		var c = this.getCanvas(event);
		
		if(c.plot.hoverIndex !== undefined){
			c.plot.draggingIndex = c.plot.hoverIndex;
			//console.log("CLICKED ON: "+c.plot.hoverIndex);
		}	
	}

	mouseUpListener(event){
		var c = this.getCanvas(event);
		
		//if(c.plot.draggingIndex !== undefined)console.log("STOPPED DRAGGING: "+c.plot.draggingIndex);
		c.plot.draggingIndex = undefined;
		
		c.plot.data.callback(c.plot.mousePoint);	
	}

	mouseMoveListener(event){
		//console.log(event);

		var c = this.getCanvas(event);
		
		if(c.plot.data.bsplineCurveFit){		
			//clear the canvas.
			c.ctx.clearRect(0, 0, c.ctx.canvas.width, c.ctx.canvas.height);
			
			//find the coordinates of the mouse pointer.
			c.plot.mousePoint.x = event.pageX-c.canvas.offsetLeft;
			c.plot.mousePoint.y = event.pageY-c.canvas.offsetTop;
			
			//does this coincide with a control point?
			c.plot.getMouseHover();
			
			if(c.plot.draggingIndex !== undefined)c.plot.dragControlPoint(c.plot.data.bsplineAllowDrag);
			
			//redraw the plot with a point at the mouse point.	
			c.plot.drawPlot();
		}	
	}

	getAxes(data){
		var xAxis = {};
		var yAxis = {};
		xAxis.min = yAxis.min = 1.E12;
		xAxis.max = yAxis.max = -1.E12;

		var j, coords, point;
		for(var i=0;i<data.series.length;i++){
			coords = data.series[i];
			for(var j=0;j<coords.length;j++){
				//allow users to specify arrays of weighted points, or just {x:, y:} objects.
				if(coords[j].p)point = coords[j].p;
				else point = coords[j];
				
				if(point.x<xAxis.min)xAxis.min = point.x;
				if(point.x>xAxis.max)xAxis.max = point.x;
				if(point.y<yAxis.min)yAxis.min = point.y;
				if(point.y>yAxis.max)yAxis.max = point.y;
			}	
		}
		
		//account for bspline control points
		if(data.bsplineCurveFit){
			for(var i=0;i<data.bsplineCurveFit.P.length;i++){
				if(data.bsplineCurveFit.P[i].p.x<xAxis.min)xAxis.min = data.bsplineCurveFit.P[i].p.x;
				if(data.bsplineCurveFit.P[i].p.x>xAxis.max)xAxis.max = data.bsplineCurveFit.P[i].p.x;
				if(data.bsplineCurveFit.P[i].p.y<yAxis.min)yAxis.min = data.bsplineCurveFit.P[i].p.y;
				if(data.bsplineCurveFit.P[i].p.y>yAxis.max)yAxis.max = data.bsplineCurveFit.P[i].p.y;
			}
		}

		if(data.reverseYAxis){
			var tmp = yAxis.max;
			yAxis.max = yAxis.min;
			yAxis.min = tmp;
		}

		return {xAxis: xAxis, yAxis: yAxis};
	}

	drawPlot(){
		var p = this.plotPanel.parameters;
		
		//draw outline around plot boundary.	
		if(this.data.hasBorder){
			this.ctx.lineWidth=1.;	
			this.ctx.strokeRect(0, 0, p.plotWidth, p.plotHeight); //black rectangular outline.
		}	
		
		//draw outline around drawing area in this plot.	
		this.ctx.lineWidth=0.5;	
		this.ctx.strokeRect(this.drawBoundary.topLeftX, this.drawBoundary.topLeftY, this.drawBoundary.width, this.drawBoundary.height); //black rectangular outline.	
		
		//draw each data series.
		if(this.data.series){
			this.drawAxisLabels(this.data, this.drawBoundary, this.ctx);
			var style;
			for(var i=0;i<this.data.series.length;i++){
				style={line:true, color:"#000000"}; //default to black line.
				if(i==1)style.color = "#0000ff"; //blue
				else if(i==2)style.color = "#00ff00";//green
				else if(i==3)style.color = "#ff0000";//red
				this.drawSeries(this.data.series[i], this.drawBoundary, this.ctx, style);
			}	
			
			//title
			if(this.data.title){
				this.ctx.textAlign = "center";
				this.ctx.font = "15px Arial";
				this.ctx.fillText(this.data.title,this.drawBoundary.topLeftX+this.drawBoundary.width/2,this.drawBoundary.topLeftY-5);
			}
		}

		if(this.data.bsplineCurveFit){
			var numPoints = 100, bsplineSeries = [], controlPointSeries=[], p;
			for(var i=0;i<numPoints;i++){
				p = this.data.bsplineCurveFit.pointOnNURBSCurve(i/(numPoints-1));
				bsplineSeries.push({x:p.x, y:p.y});
			}
			//draw the bspline curve fit with a red line.
			this.drawSeries(bsplineSeries, this.drawBoundary, this.ctx, {line:true, filledPoint: false, color:'#ff0000'});
			//draw control points as green boxes, control polygon with green line.
			for(var i=0;i<this.data.bsplineCurveFit.P.length;i++)controlPointSeries.push(this.data.bsplineCurveFit.P[i].p);
			this.drawSeries(controlPointSeries, this.drawBoundary, this.ctx, {line:true, filledPoint: true, color:'#00ff00'});		
			//if the mouse is hovering over a control point, outline it in yellow.
			if(this.hoverIndex !== undefined){
				var style = {filledPoint: false, color:'#ffff00', width:15, lineWidth:2};
				this.drawMousePoint(this.screenControlPoints[this.hoverIndex], this.ctx, style);
			}
			/*
			else if(this.mousePoint){
				//draw a yellow box that follows the mouse point.
				this.drawMousePoint(this.mousePoint, this.ctx, {filledPoint: true, color:'#777700'});
			}
			*/
		}

	}

	drawSeries(series, drawBoundary, ctx, style){
		if(!style)style = {line: true, color:'#0000000'};

		ctx.save();
		ctx.strokeStyle=style.color;
		ctx.fillStyle=style.color;	
		ctx.lineWidth=1.5;
		ctx.beginPath();

		var coord;
		for(var i=0; i<series.length;i++){
			//allow user to specify series as array of weighted points, or just {x:, y:} objects.
			if(series[i].p)coord = this.toScreenCoord(series[i].p)
			else coord = this.toScreenCoord(series[i]);
			
			if(style.filledPoint)ctx.fillRect(coord.x-2.5, coord.y-2.5, 5, 5);
			if(style.line)ctx.lineTo(coord.x, coord.y);
		}	
		ctx.stroke();
		ctx.closePath();
		ctx.restore();
	}

	drawMousePoint(point, ctx, style){
		if(!style)style = {color:'#0000000'};

		var boxWidth = 5;
		if(style.width)boxWidth = style.width;

		ctx.save();	
		ctx.strokeStyle=style.color;
		ctx.fillStyle=style.color;
		if(style.lineWidth)ctx.lineWidth = style.lineWidth;
		if(style.filledPoint)ctx.fillRect(point.x-boxWidth/2, point.y-boxWidth/2, boxWidth, boxWidth);
		else ctx.strokeRect(point.x-boxWidth/2, point.y-boxWidth/2, boxWidth, boxWidth);
		ctx.restore();
	}

	drawAxisLabels(data, drawBoundary, ctx){
		ctx.lineWidth=1.5;
		//tics
		var numXdivisions = 10;
		var numYdivisions = 10;
		var numXtics = numXdivisions + 1;
		var numYtics = numYdivisions + 1;		
		if(data.numXtics)numXtics = data.numXtics;
		else if(data.xTicStep)numXtics = ((drawBoundary.xAxis.max - drawBoundary.xAxis.min)/data.xTicStep).round(0) + 1;
		if(data.numYtics)numYtics = data.numYtics;
		else if(data.yTicStep)numYtics = ((drawBoundary.yAxis.max - drawBoundary.yAxis.min)/data.yTicStep).round(0) + 1;
		var step_xLabel = (drawBoundary.xAxis.max - drawBoundary.xAxis.min)/(numXtics-1);
		var step_yLabel = (drawBoundary.yAxis.max - drawBoundary.yAxis.min)/(numYtics-1);		
		var ticLength = drawBoundary.width/50;
		//
		ctx.font = "12px Arial";
		//
		var x, y, label;
		
		//X axis
		var dec = Math.log10(step_xLabel);
		if(dec>=0)dec = 1;
		else dec = Math.abs(dec).round(0);
		ctx.beginPath();
		for(var i=0;i<numXtics;i++){
			label = (drawBoundary.xAxis.min+(step_xLabel*i)).round(dec);		
			x = this.toScreenCoord({x:label}).x;
			y = drawBoundary.topLeftY + drawBoundary.height;
			ctx.moveTo(x,y);
			
			//draw tic label
			ctx.textAlign = "center";
			ctx.fillText(label,x,y+15);
			
			//draw tic
			y -= ticLength;			
			ctx.lineTo(x,y);
		}
		ctx.stroke();
		ctx.closePath();
		if(data.xAxisLabel){
			ctx.textAlign = "center";
			ctx.font = "13px Arial";
			x = drawBoundary.topLeftX + drawBoundary.width/2;
			y = drawBoundary.topLeftY + drawBoundary.height + 35;
			ctx.fillText(data.xAxisLabel, x, y);
		}
		
		//Y axis
		dec = Math.log10(Math.abs(step_yLabel));
		if(dec>=0)dec = 1;
		else dec = Math.abs(dec).round(0);
		ctx.beginPath();
		for(var i=0;i<numYtics;i++){
			x = drawBoundary.topLeftX;
			label = (drawBoundary.yAxis.min+(step_yLabel*i)).round(dec);
			y = this.toScreenCoord({y:label}).y;
			ctx.moveTo(x,y);
			
			//draw tic label
			ctx.textAlign = "right";
			ctx.fillText(label,x-5,y);
			
			//draw tic			
			x += ticLength;			
			ctx.lineTo(x,y);
		}
		ctx.stroke();
		ctx.closePath();
		if(data.yAxisLabel){
			ctx.font = "13px Arial";			
			x = drawBoundary.topLeftX - 35;
			y = drawBoundary.topLeftY + drawBoundary.height/2;
 
			ctx.save();
			ctx.translate(x, y);
			ctx.rotate(-Math.PI/2);
			ctx.textAlign = "center";
			ctx.fillText(data.yAxisLabel, 0, 0);
			ctx.restore();		
		}
	}
}

/////
/////
/*
PlotPanel.prototype.makeTestPlot = function(){
	var parameters = {
		plotContainer: document.getElementById('plotContainer'),
		plotWidth: 350,
		plotMargin: {left:50, right:35, top:35, bottom:50},
		numPlotsAcrossPage: 2,
		hasBorder: false,
	};	
	
	var plotPanel = new PlotPanel(parameters);
	
	//
	var dataSeries = function(numPoints, theFunction){
		var series = new Array(numPoints);

		var range = 720.;
		var step = range/numPoints;			
		
		var x, y;
		for(var i=0;i<=numPoints;i++){
			x = -(range/2.) + i*step;
			y = theFunction(x);
			series[i] = {x:x, y:y};
		}

		return series;
	}
	
	//
	var cosPlot = {
		coords:dataSeries(100, function(x){return Math.cos(x*Math.PI/180.)}), 
		xTicStep:90, 
		title:"cosine", 
		xAxisLabel:"\u03b8 (degrees)", 
		yAxisLabel:"Cosine(\u03b8)",
		//hasBorder: true,
	};
	plotPanel.addPlot(cosPlot);


	//
	var sinPlot = {
		coords:dataSeries(100, function(x){return Math.sin(x*Math.PI/180.)}), 
		xTicStep:90, 
		title:"sine", 
		xAxisLabel:"\u03b8 (degrees)", 
		yAxisLabel:"Sine(\u03b8)",
		//hasBorder: true,		
	};
	plotPanel.addPlot(sinPlot);
	//
	var tanPlot = {
		coords:dataSeries(100, function(x){return Math.tan(x*Math.PI/180.)}), 
		xTicStep:90, 
		title:"tan", 
		xAxisLabel:"\u03b8 (degrees)", 
		yAxisLabel:"Tan(\u03b8)",
		//hasBorder: true,
	};
	plotPanel.addPlot(tanPlot);
	//
	var sqrtPlot = {
		coords:dataSeries(100, Math.sqrt), 
		xTicStep:90, 
		title:"sqrt", 
		xAxisLabel:"\u03b8 (degrees)", 
		yAxisLabel:"\u221a\u03b8",
		//hasBorder: true,
	};
	plotPanel.addPlot(sqrtPlot);
	//
	
	var sketchTextDiv = document.getElementById('bladeSketchTextContainer');
	sketchTextDiv.innerHTML = "<p>BLAH!</p>"
	//console.log(sketchTextDiv.innerHTML);
}
*/
//
Number.prototype.round = function(places) {
  var result = +(Math.round(this + "e+" + places)  + "e-" + places);
  if(isNaN(result))return 0.;
  return result;
}	