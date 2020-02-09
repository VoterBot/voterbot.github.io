"use strict";
class DrawModel{
	constructor(canvas_width, canvas_height, canvasDiv, isDesignMode){
		this.canvas_width = canvas_width;
		this.canvas_height = canvas_height;
	
		// Three.js rendering basics.
		this.scene;
		this.camera; 
		this.renderer;  
		this.controls;
	
		this.cameraAndLight;  // Object holding both camera and light, for rotation.

		//////////////////
		//initialise.
		try{
			this.renderer = new THREE.WebGLRenderer({antialias:true});
			//this.renderer.setPixelRatio(window.devicePixelRatio*2); //fixes an antialiasing problem in firefox.			
			//this.renderer.setPixelRatio(1.25); //fixes an antialiasing problem in firefox.			
			this.renderer.setSize(canvas_width, canvas_height);
			canvasDiv.appendChild(this.renderer.domElement);
			
			if (!this.renderer) {
				canvasDiv.innerHTML = "Sorry, WebGL is required but is not available.";
				return;
			}
			var background = 0x444444;  // dark grey background
			if(isDesignMode) background = 0xffffff;  // white background
			
			this.renderer.setClearColor(background);
			
			this.camera = new THREE.PerspectiveCamera(20, canvas_width/canvas_height, 0.1, 1000);
			
			//initialise the scene.
			this.scene = new THREE.Scene();		
			this.createWorld(isDesignMode);
			
			//add mouse controls for rotating, panning and zooming
			this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement );
			this.controls.addEventListener('change', this.render);
			//this.controls.damping = 0.2;
			
			//point to the the current parent object in the OrbitControls object.  render() needs it.
			this.controls.drawModelObject = this;
		}
		catch(e){
			document.write("Sorry, an error occurred: " + e);
			return;
		}
	}
	
	//utility methods.
	getByName(name){
		var results = [];
		for(var i=0;i<this.scene.children.length;i++){
			if(this.scene.children[i].name.indexOf(name)>=0)results.push(this.scene.children[i]);
		}	
		
		return results;
	}
	
	removeFromScene(name){
		this.scene.remove(this.scene.getObjectByName(name));
	}
	
	render(){
		//when render() is used as a callback to OrbitControls, make 'this' refer to calling DrawModel object.
		var drawModelObject = this;
		if(this.drawModelObject)drawModelObject = this.drawModelObject;
		drawModelObject.renderer.render(drawModelObject.scene, drawModelObject.camera);
	}
	
	//methods
	createWorld(isDesignMode){
		/* Add the camera and a light to the scene, linked into one object. */
		var light1 = new THREE.DirectionalLight();
		light1.position.set(0,500,500);
		var light2 = new THREE.DirectionalLight();
		light2.position.set(0,-500,-500);
		var light3 = new THREE.AmbientLight( 0x202020 ); // soft white light						

		this.cameraAndLight = new THREE.Object3D();	
		this.cameraAndLight.add(this.camera);	
		this.cameraAndLight.add(light1);
		this.cameraAndLight.add(light2);		
		this.cameraAndLight.add(light3);

		//default position
		var position = {x:-12, y:-0.77, z:12};
		if(isDesignMode){
			position = {x: 0., y: 0., z:2.};

			//translate and rotate the cameraAndLight object
			this.cameraAndLight.translateY(1.25);
			this.cameraAndLight.rotateOnAxis(new THREE.Vector3( 0, 0, 1 ), Math.PI/2);
			this.cameraAndLight.rotateOnAxis(new THREE.Vector3( 1, 0, 0 ), -Math.PI/3);		
		}	
		
		this.camera.position.set(position.x, position.y, position.z);
		//this.camera.rotation.set(x, y, z);
		
		//add cameraAndLight object to scene.
		this.scene.add(this.cameraAndLight);		
	}

	addModelToScene(model, displayHasChanged){
		var numUSteps = 60;
		var numVSteps = 40;		
		
		//input: model -> SolidObject; 
		//input: displayHasChanged: (optional) -> "Display" menu item picked, so only update that item here.
		
		////
		//outer NURBS surface.
		var params = {
			model:model.s, 
			name: "model.outerSurface",
			isVisible: model.showOuterSurface,
			hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showOuterSurface"),
			surfaceIndex:0, 
			numUSteps:numUSteps, 
			numVSteps:numVSteps,
		};
		if(this.toggleDisplay(params, "drawSurfaceMesh"))return;	

		//inner NURBS surfaces.
		if(model.s.curves.length>1){
			params = {
				model:model.s,
				name: "model.innerSurface",
				isVisible: model.showInnerSurface, 
				hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showInnerSurface"),
				numUSteps:numUSteps, 
				numVSteps:numVSteps,
			};

			if(this.toggleDisplay(params, "drawInnerSurfaces"))return;
		}

		////
		//surface derivatives.
		params = {
			model:model.s, 
			name: "model.surfaceDerivatives", 
			isVisible: model.showSurfaceDerivatives, 
			hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showSurfaceDerivatives"),
			numUSteps:20, 
			numVSteps:20,
		};
		if(this.toggleDisplay(params, "drawSurfaceDerivatives"))return;
		
		////
		//surface panels for the panel solver.
		if(model.s.mesh.bladePanel){
			params = {
				mesh: model.s.mesh,
				isBladeMesh: true,
				name: "model.bladeMesh",
				isVisible: model.showBladePanelMesh,
				hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showBladePanelMesh"),
			};
			if(this.toggleDisplay(params, "drawPanelMesh"))return;
			
			//coordinate systems.
			params.isVisible = model.showBladePanelCoordinateSystems;
			params.name = "model.bladeCS";
			params.hasChanged = (displayHasChanged === undefined)?undefined:(displayHasChanged === "showBladePanelCoordinateSystems");
			if(this.toggleDisplay(params, "drawPanelCoordinateSystems"))return;		
		}
		//wake panels.
		if(model.s.mesh.wakePanel){
			params = {
				mesh: model.s.mesh,
				isWakeMesh: true,
				name: "model.wakeMesh",
				isVisible: model.showWakePanelMesh,
				hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showWakePanelMesh"),
			};
			if(this.toggleDisplay(params, "drawPanelMesh"))return;
			
			//coordinate systems.
			params.isVisible = model.showWakePanelCoordinateSystems;
			params.name = "model.wakeCS";
			params.hasChanged = (displayHasChanged === undefined)?undefined:(displayHasChanged === "showWakePanelCoordinateSystems");
			if(this.toggleDisplay(params, "drawPanelCoordinateSystems"))return;
		}
		
		////
		//curves
		params = {
			model: model.s,
			name: "model.curves",
			isVisible: model.showAllCurves,
			hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showAllCurves"),
		};
		if(this.toggleDisplay(params, "drawAllCurves"))return;	
		
		//axes.
		params = {
			model: model,
			name: "model.axes",
			isVisible: model.showAxes,
			hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showAxes"),		
		};
		if(this.toggleDisplay(params, "drawAxes"))return;	
		
		//and draw everything.
		this.render();		
	}

	toggleDisplay(p, draw){
		//if an item on the "Display" menu has been selected, toggle visibility for all items with that name.
		if(p.hasChanged){
			var x = this.getByName(p.name);
			if(x && x.length>0){
				for(var i=0;i<x.length;i++)
					x[i].visible = p.isVisible;
			}	
			else p.isNew = true; //threejs object doesn't exist yet.  Flag it for redraw.
		}	
		
		//if model has changed, or if the threejs object for it hasn't been created yet, create it now.
		if(p.hasChanged === undefined || p.isNew)this[draw](p);	

		//stop updating things if this is a "Display" menu pick.
		if(p.hasChanged){
			this.render();
			return true;
		}
		return false;
	}

	addBladeDesignToScene(model){
		//outer NURBS surface.
		var params = {
			model:model.s, 
			isVisible: true, 
			name: "designSurface.surface", 
			surfaceIndex:0, 
			numUSteps:60, 
			numVSteps:40,
			surfaceColor:0x000001,
			opacity:0.4,
		};
		this.drawSurfaceMesh(params);
		
		////
		//curves
		params = {
			model: model.s,
			name: "designSurface.curves",
			isVisible: true,
		}
		this.drawAllCurves(params);
		
		//axes.
		params = {
			model: model,
			name: "designSurface.axes",
			isVisible: true,
			hideAxisPlanes: true,
		}
		this.drawAxes(params);

		//and draw everything.
		this.render();		
	}

	addSolverResultsToScene(model, gradientCanvas){
		this.createScalarFieldShading(model, gradientCanvas);	

		//surface panels for the panel solver.
		var params;
		if(model.s.mesh.bladePanel){
			params = {
				mesh: model.s.mesh,
				isResults: true,
				isScalarField: true,
				name: "results.bladeMesh",
				isVisible: true,
			}
			this.drawPanelMesh(params);
		}
		
		//axes.
		params = {
			model: model,
			name: "results.axes",
			isVisible: true,
			hideAxisPlanes: true,
		}
		this.drawAxes(params);
		
		//and draw everything.
		this.render();		
	}

	getScalarFieldColourMap(){
		//linear colour map between 100% and 0%.  Maps to HSL colour values 0 (red) to blue (240).
		var colourMap = {
			stop: undefined,
			maxHue: 0,
			minHue: 0,
			
			set: function(){
				//non-linear map to de-emphasise green.
				var x = new Array(5);
				x[0] = {val: 1., hue: 0};//red.
				x[1] = {val: 0.666, hue: 60};//yellow.
				x[2] = {val: 0.5, hue: 120};//green.
				x[3] = {val: 0.333, hue: 180};//cyan.
				x[4] = {val: 0., hue: 240};//blue.

				this.stop = x;
				
				this.maxHue = x[x.length-1].hue;
				this.minHue = x[0].hue;
				
				for(var i=0; i<x.length; i++)x[i].hslString = this.getHSLString(x[i].val);			
			},
			
			getHue: function(x){
				if(!this.stop)this.set();
				x = x.round(6);
				//x: [1., 0.]
				//linear interpolation between the colour stops.
				var percent;
				for(var i=1;i<this.stop.length;i++){
					if(x<=this.stop[i-1].val && x>=this.stop[i].val){
						percent = (x-this.stop[i].val)/(this.stop[i-1].val - this.stop[i].val);
						return (this.stop[i].hue - percent*(this.stop[i].hue - this.stop[i-1].hue)).round(0);
					}
				}
				
				return this.maxHue;
			},
			
			getHSLString: function(x){
				return ("hsl("+this.getHue(x)+", 100%, 50%)");
			},
		};
		
		colourMap.set();
		
		return colourMap;
	}

	createScalarLegendBar(colourMap, min, max, gradientCanvas, units){
		///////
		//http://www.w3schools.com/tags/canvas_createlineargradient.asp
		var legendWidth = 100;
		var bar = {
			width: 20,
			height: 200,
			left: 80,
			top: 15,
		};

		var g = gradientCanvas.style;
		g.display = "block";
		g.position = "absolute";	
		g.zIndex = 2;
		g.left = (this.canvas_width - legendWidth)+"px";
		g.top = "40px";
		//
		gradientCanvas.height = bar.height+bar.top; //ensure canvas is large enough.
		
		var ctx=gradientCanvas.getContext("2d");
		
		//create a linear gradient along a vertical line, lower values at the bottom.
		var grd=ctx.createLinearGradient(bar.left,bar.top+bar.height,bar.left,bar.top);

		//make colour map.
		for(var i=0;i<colourMap.stop.length;i++){
			grd.addColorStop(colourMap.stop[i].val, colourMap.stop[i].hslString);	
		}
		
		ctx.fillStyle=grd;
		ctx.fillRect(bar.left, bar.top, bar.width, bar.height);	
		/*
		ctx.strokeStyle = "black";
		ctx.lineWidth = 0.2;
		ctx.strokeRect(bar.left, 0, bar.width, 200);
		*/
		
		///////////////////////
		///////////////////////
		//draw max, min and zero tics.
		ctx.strokeStyle="black";
		ctx.fillStyle="black";	
		ctx.font="10px Arial";
		
		//max
		var x = bar.left+bar.width+5;
		var y = bar.top+10;
		
		var label = max.round(2);	
		ctx.fillText(label,x,y);
		///////
		///////
		
		//min
		x = bar.left+bar.width+5;
		y = bar.top+bar.height;
		
		label = min.round(2);	
		ctx.fillText(label,x,y);
		///////
		///////
		
		//zero
		if(max >= 0. && min <= 0.){
			ctx.lineWidth = 0.5;
			
			var p = -min/(max-min);
			var colour = colourMap.getHSLString(p);
			
			ctx.strokeStyle=colour;
			ctx.fillStyle=colour;	
			
			y = bar.top+(bar.height*(1.-p));
			
			x = bar.left;
			ctx.moveTo(x,y);	
			x -= 5
			ctx.lineTo(x,y);

			x = bar.left+bar.width+5;
			y += 5;
			
			label = "0.0";	
			ctx.fillText(label,x,y);
			ctx.stroke();
			ctx.closePath();
		}	
		///////
		///////

		//units.
		if(units){
			ctx.fillStyle="black";	
			ctx.textAlign = "center";		
			
			y = bar.top - 5;
			x = bar.left+bar.width/2;
			ctx.fillText(units,x,y);
		}	

	}

	createScalarFieldShading(model, gradientCanvas){
		//shade each blade panel node with a colour based on field point value at collocation point.

		if(!model.showResults)model.set_Cp(); //default to a display of pressure coefficient.
		var scalarField = model.showResults.field;

		//colour map.
		var colourMap = this.getScalarFieldColourMap();
		
		var nPanels = model.s.mesh.nPanels, mPanels = model.s.mesh.mPanels;
		var numPanels = nPanels*mPanels, panels = model.s.mesh.bladePanel;
		
		//deal with scalar field values on the key blade.
		var max = -1.E6, min= 1.E6; 
		for(var i=0;i<numPanels;i++){
			//if(panels[i].r_on_R > 0.15){
				if(panels[i][scalarField] > max)max = panels[i][scalarField];
				if(panels[i][scalarField] < min)min = panels[i][scalarField];
			//}	
		}
		var range = max - min;
		
		//create a coloured legened bar depicting this range.
		this.createScalarLegendBar(colourMap, min, max, gradientCanvas, model.showResults.units);

		//average collocation field point at adjoining nodes; use this value to colour each node.
		var a,b,c,d, m, n, index=-1, val, colour;
		for(m=0;m<mPanels;m++){	
			for(n=0;n<nPanels;n++){
				index++;
		
				//wash out bad values.
				if(panels[index][scalarField] > max || panels[index][scalarField] < min)panels[index].luminance = "100%";
				else panels[index].luminance = "50%";	
				
				a = index;
				b = index+1; 
				c = b-nPanels;
				d = c-1;
				
				//first blade element.
				if(m==0){
					if(n<nPanels-1){
						val = (panels[a][scalarField] + panels[b][scalarField])/2.;
						colour = colourMap.getHue((val - min)/range);
						panels[a].node[0].colour = colour;
						panels[b].node[3].colour = colour;
					}
					
					//trailing edge.
					if(n==0)panels[a].node[3].colour = colourMap.getHue((panels[a][scalarField] - min)/range);
					if(n==nPanels-1)panels[a].node[0].colour = colourMap.getHue((panels[a][scalarField] - min)/range);			
				}
				else{
					if(n<nPanels-1){
						val = (panels[a][scalarField] + panels[b][[scalarField]] + panels[c][scalarField] + panels[d][scalarField])/4.;
						colour = colourMap.getHue((val - min)/range);
						
						panels[a].node[0].colour = colour;
						panels[b].node[3].colour = colour;
						panels[c].node[2].colour = colour;
						panels[d].node[1].colour = colour;
					}
					
					//trailing edge.
					if(n==0){
						val = (panels[a][scalarField] + panels[d][scalarField])/2.;
						colour = colourMap.getHue((val - min)/range);
					
						panels[a].node[3].colour = colour;
						panels[d].node[2].colour = colour;				
					}	
					if(n==nPanels-1){
						val = (panels[a][scalarField] + panels[d][scalarField])/2.;
						colour = colourMap.getHue((val - min)/range);
					
						panels[a].node[0].colour = colour;
						panels[d].node[1].colour = colour;				
					}	

				}

				//last row.
				if(m==mPanels-1){
					if(n<nPanels-1){
						val = (panels[a][scalarField] + panels[b][scalarField])/2.;
						colour = colourMap.getHue((val - min)/range);
						panels[a].node[1].colour = colour;
						panels[b].node[2].colour = colour;
					}
					
					//trailing edge.
					if(n==0)panels[a].node[2].colour = colourMap.getHue((panels[a][scalarField] - min)/range);
					if(n==nPanels-1)panels[a].node[1].colour = colourMap.getHue((panels[a][scalarField] - min)/range);
				}
			}	
		}
	}

	drawInnerSurfaces(p){
		/*
		params = {
			model:model.s,
			name: "model.innerSurface",
			isVisible: model.showInnerSurface, 
			hasChanged: (displayHasChanged === undefined)?undefined:(displayHasChanged === "showInnerSurface"),
			numUSteps:numUSteps, 
			numVSteps:numVSteps,
		};
		*/
		
		for(var i=1;i<p.model.curves.length;i++){
			p.name = "model.innerSurface."+i;
			p.surfaceIndex = i;
			this.drawSurfaceMesh(p);
		} 
		
		//a centreline for the final, innermost curve.
		p.name = "model.innerSurface.centreline";
		this.drawCentreLine(p);
	}

	drawAllCurves(p){
		/*
		params = {
			model: model.s,
			name: "model.curves",
			isVisible: model.showAllCurves,
		}
		*/	

		this.removeFromScene(p.name);
		if(!p.isVisible)return;

		var curveContainer = new THREE.Object3D();
		curveContainer.name = p.name;
		
		var surface = p.model;
		//use the chord length of the first curve to scale the size of control point spheres
		var sphereRadius = surface.curves[0][0].scale/75;

		//draw individual surface curves
		var j, params;
		for(var i=0;i<surface.curves[0].length;i++){
			//add control points and polygon
			params = {
				curve: surface.curves[0][i].c,			
				sphereRadius: sphereRadius,
				isVisible: (p.isVisible && surface.curves[0][i].showControlPoints),
				name: p.name+".surfaceCurve."+i,
			};
			curveContainer.add(this.getControlPoints(params));

			//add the curve itself
			params.isVisible = (p.isVisible && surface.curves[0][i].showCurve);
			params.isSurfaceCurve = true;
			curveContainer.add(this.getCurve(params));
			//add any inner geometry curves
			params.isSurfaceCurve = false;		
			for(j=1;j<surface.curves.length;j++){
				params.name = p.name+".innerCurve."+i+"."+j;
				params.curve = surface.curves[j][i].c;
				curveContainer.add(this.getCurve(params));
			}	
		}
		
		this.scene.add(curveContainer);
	}

	getControlPoints(p){
		/*
		params = {
			curve: surface.curves[0][i].c,			
			sphereRadius: sphereRadius,
			isVisible: (p.isVisible && surface.curves[0][i].showControlPoints),
			name: p.name+".surfaceCurve."+i,
		};
		*/
		
		var curve = p.curve;

		var name_spheres = p.name+".controlPoints";
		var name_line = p.name+".controlPolygon";
		var curveContainer = new THREE.Object3D();
		
		var spheres = new THREE.Object3D();
		spheres.name = name_spheres;
		spheres.visible = p.isVisible;

		//add control points as little red spheres
		var geom = new THREE.SphereGeometry(p.sphereRadius,6,6); //radius, width segments, height segments
		var material = new THREE.MeshPhongMaterial({color: 0xff0000, specular:0x333333, shininess: 100});

		while(spheres.children.length < curve.P.length)spheres.children.push(new THREE.Mesh(geom, material));
		
		material.dispose();
		geom.dispose();

		for (var i=0; i<curve.P.length; i++){
			spheres.children[i].position.set(curve.P[i].p.x, curve.P[i].p.y, curve.P[i].p.z);					
		}
		curveContainer.add(spheres);	

		////
		//add control polygon as blue lines
		geom = new THREE.Geometry();
		
		while(geom.vertices.length < curve.P.length)geom.vertices.push(new THREE.Vector3(0.,0.,0.));
		
		for (var i=0; i<curve.P.length; i++){
			geom.vertices[i].set(curve.P[i].p.x, curve.P[i].p.y, curve.P[i].p.z);
		}	
		
		var material = new THREE.LineBasicMaterial({color: 0x0000ff});	
		var line = new THREE.Line(geom, material);
		line.name = name_line;
		line.visible = p.isVisible;		

		material.dispose();
		geom.dispose();

		curveContainer.add(line);
		
		return curveContainer;
	}

	getCurve(p){
		/*
		params = {
			curve: surface.curves[0][i].c,			
			isVisible: (p.isVisible && surface.curves[0][i].showCurve),
			name: p.name+".surfaceCurve."+i,
			isSurfaceCurve = true,
		};
		*/
		
		var curve = p.curve;
		
		//add the curve itself as a green line
		var geom = new THREE.Geometry();
		
		var numSegments = 100;
		
		var uStep = (1./numSegments).round(8);	
		var z = curve.P[0].p.z, u, point;		
		for(var i=0; i<=numSegments;i++){
			u = (i*uStep).round(8);
			if(u>1.)u=1.;
			point = curve.pointOnNURBSCurve(u);
			if(point.z==undefined)geom.vertices.push(new THREE.Vector3(point.x, point.y, z));
			else geom.vertices.push(new THREE.Vector3(point.x, point.y, point.z));
		}	
		
		var curveColour = 0x00ff00; //green for surface curves
		if(!p.isSurfaceCurve)curveColour = 0xff0000; //red for inner curves
		var material = new THREE.LineBasicMaterial({color: curveColour});
		
		var line = new THREE.Line(geom, material);
		line.name = p.name;
		line.visible = p.isVisible;		
		
		material.dispose();
		geom.dispose();

		return line;
	}

	drawSurfaceDerivatives(p){	
		/*
		params = {
			model:model.s, 
			isVisible: model.showSurfaceDerivatives, 
			name: "model.surfaceDerivatives", 
			numUSteps:20, 
			numVSteps:20,	
		};
		*/

		var surface = p.model;
		
		var name_u = p.name+".u";
		var name_v = p.name+".v";
		var name_n = p.name+".n";
		
		this.removeFromScene(name_u);
		this.removeFromScene(name_v);
		this.removeFromScene(name_n);
		if(!p.isVisible)return;
		
		//one geometry for each local axis, containing all vectors in that direction.
		
		var geom_u = new THREE.Geometry();
		var geom_v = new THREE.Geometry();	
		var geom_n = new THREE.Geometry();	
		
		//if the original curves on which the solid is constructed are defined with anti-clockwise
		//numbered control points, the outward normal is du.cross.dv.  If clockwise, n = dv.cross.du.
		var hasReversedNormal = false;
		if(surface.outwardNormalDirection==-1)hasReversedNormal = true;
		
		var uStep = (1./p.numUSteps);
		var vStep = (1./p.numVSteps);
		
		var i, j, d=1, u, v, derivs, p0, p_du, p_dv, p_n, e, f, g, h;
		for(i=0;i<=p.numUSteps;i++){
			u = (i*uStep); if(u>1.)u=1.;
			for(j=0;j<=p.numVSteps;j++){
				v = (j*vStep); if(v>1.)v=1.;
				
				derivs = surface.rationalSolidDerivs(u, v, 0, d);
				p0 = derivs[0][0][0].p;
				p_du = derivs[1][0][0].p;
				p_dv = derivs[0][1][0].p;

				//calculate surface normal
				e = p_du.dot(p_du);
				f = p_dv.dot(p_dv);
				g = p_du.dot(p_dv);
				h = Math.sqrt(e*f - g*g);
				if(hasReversedNormal)p_n = p_dv.cross(p_du).divide(h);
				else p_n = p_du.cross(p_dv).divide(h);	//default.			

				//scale these vectors so they look nice on the screen
				p_du = p_du.divide(25);
				p_dv = p_dv.divide(25);
				//p_dw = p_dw.divide(25);				
				p_n = p_n.divide(25);				
				p_du = p_du.plus(p0);
				p_dv = p_dv.plus(p0);
				//p_dw = p_dw.plus(p0);				
				p_n = p_n.plus(p0);
				
				//u direction
				geom_u.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
				geom_u.vertices.push(new THREE.Vector3(p_du.x, p_du.y, p_du.z));
				//v direction
				geom_v.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
				geom_v.vertices.push(new THREE.Vector3(p_dv.x, p_dv.y, p_dv.z));
				//normal
				geom_n.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
				geom_n.vertices.push(new THREE.Vector3(p_n.x, p_n.y, p_n.z));
			}
		}     

		var material_u = new THREE.LineBasicMaterial({color: 0x0000ff}); //blue for u direction
		var material_v = new THREE.LineBasicMaterial({color: 0x00ff00}); //green for v direction
		var material_n = new THREE.LineBasicMaterial({color: 0xff0000}); //red for surface normal
		
		//create Three.LineSegments objects, each containing all vectors in one parametric direction.
		var line_u = new THREE.LineSegments(geom_u, material_u);
		var line_v = new THREE.LineSegments(geom_v, material_v);
		var line_n = new THREE.LineSegments(geom_n, material_n);
		
		line_u.name = name_u;
		line_v.name = name_v;
		line_n.name = name_n;
		
		line_u.visible = line_v.visible = line_n.visible = p.isVisible;		
		
		this.scene.add(line_u);
		this.scene.add(line_v);
		this.scene.add(line_n);
		
		geom_u.dispose();
		geom_v.dispose();
		geom_n.dispose();
		material_u.dispose();
		material_v.dispose();
		material_n.dispose();
	}

	drawSurfaceMesh(p){	
		/*
		var params = {
			model:model.s, 
			isVisible: model.showOuterSurface, 
			name: "model.outerSurface", 
			surfaceIndex:0, 
			numUSteps:numUSteps, 
			numVSteps:numVSteps,	
		};
		*/
		
		this.removeFromScene(p.name);
		if(!p.isVisible)return;

		var surface = p.model;

		var surfaceGeometry = new THREE.Geometry();
		
		var uStep = (1./p.numUSteps); //javascript floating point precision needs to be tidied up.
		var vStep = (1./p.numVSteps);		

		var i, j, u, v, p0, p1, p2, p3, k=0;
		for(i=0;i<p.numUSteps;i++){
			u = i*uStep; if(u>1.)u=1.;
			for(j=0;j<p.numVSteps;j++){
				v = (j*vStep); if(v>1.)v=1.;
				
				p0 = surface.surfacePoint(u, v, p.surfaceIndex);
				p1 = surface.surfacePoint(u+uStep, v, p.surfaceIndex);
				p2 = surface.surfacePoint(u+uStep, v+vStep, p.surfaceIndex);
				p3 = surface.surfacePoint(u, v+vStep, p.surfaceIndex);
				
				//quad vertices
				surfaceGeometry.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z)); //0
				surfaceGeometry.vertices.push(new THREE.Vector3(p1.x, p1.y, p1.z)); //1
				surfaceGeometry.vertices.push(new THREE.Vector3(p2.x, p2.y, p2.z)); //2
				surfaceGeometry.vertices.push(new THREE.Vector3(p3.x, p3.y, p3.z)); //3
				
				//two triangular faces
				surfaceGeometry.faces.push(new THREE.Face3(k, k+1, k+3));
				surfaceGeometry.faces.push(new THREE.Face3(k+3, k+1, k+2));			
				//
				
				k+=4;	
			}
		}     

		surfaceGeometry.computeFaceNormals();
		
		if(p.surfaceIndex>0)p.surfaceColor = 0xff0000; //red for interior surfaces
		if(!p.surfaceColor)p.surfaceColor = 0xffffff; //white	
		
		var material =	new THREE.MeshPhongMaterial( {
			color: p.surfaceColor,
			specular:0x333333,
			shininess: 100,
			side: THREE.DoubleSide,
		});
		if(p.opacity){
			material.transparent = true;
			material.opacity = p.opacity;
		}

		var surfaceMesh = new THREE.Mesh(surfaceGeometry, material);
		surfaceMesh.name = p.name;
		surfaceMesh.visible = p.visible;		
		this.scene.add(surfaceMesh);

		surfaceGeometry.dispose();
		material.dispose();
	}
	
	drawPanelMesh(p){	
		/*
		params = {
			mesh: model.s.mesh,
			isBladeMesh: true,
			name: "model.bladeMesh",
			isVisible: model.showBladePanelMesh,
		}
		*/

		var name_surfaceMesh = p.name+".surfaceMesh";
		var name_meshLines = p.name+".meshLines";
		
		this.removeFromScene(name_surfaceMesh);
		this.removeFromScene(name_meshLines);	

		if(!p.isVisible)return; //no need to update geometry if it's hidden.	

		var surfaceColor, transparent, opacity, panelMesh, numPanels;
		if(p.isBladeMesh || p.isResults){
			panelMesh = p.mesh.bladePanel;
			surfaceColor = 0xcae1ff; //light steelblue 1; for blade panels.
			transparent = false;
			opacity = 1.;
			numPanels = p.mesh.nPanels * p.mesh.mPanels * p.mesh.numBlades;
			if(p.isResults)numPanels = p.mesh.nPanels * p.mesh.mPanels; //single blade for results display.
		}	
		else if(p.isWakeMesh){
			panelMesh = p.mesh.wakePanel;		
			surfaceColor = 0xffff00; //yellow; for wake panels.
			transparent = true;
			opacity = 0.3;
			numPanels = p.mesh.mPanels * p.mesh.wPanels * p.mesh.turns * p.mesh.numBlades;		
		}

		var surfaceGeometry = new THREE.Geometry();
		var outlineGeom = new THREE.Geometry();

		var j, k=0, v, faceA, faceB, luminance, hsl = new Array(4);
		for(var i=0;i<numPanels;i++){
			//quad vertices
			for(j=0;j<4;j++){
				surfaceGeometry.vertices.push(new THREE.Vector3(panelMesh[i].node[j].x, panelMesh[i].node[j].y, panelMesh[i].node[j].z));
			}

			//two triangular faces
			faceA = new THREE.Face3(k, k+3, k+1);
			faceB = new THREE.Face3(k+1, k+3, k+2);
			
			//colour each vertex.
			if(p.isScalarField){
				luminance = panelMesh[i].luminance;
				
				hsl[0] = new THREE.Color("hsl("+panelMesh[i].node[0].colour+", 100%, "+luminance+")");
				hsl[1] = new THREE.Color("hsl("+panelMesh[i].node[1].colour+", 100%, "+luminance+")");
				hsl[2] = new THREE.Color("hsl("+panelMesh[i].node[2].colour+", 100%, "+luminance+")");
				hsl[3] = new THREE.Color("hsl("+panelMesh[i].node[3].colour+", 100%, "+luminance+")");			
					
				faceA.vertexColors[0] = hsl[0];
				faceA.vertexColors[1] = hsl[3];
				faceA.vertexColors[2] = hsl[1];

				//
				faceB.vertexColors[0] = hsl[1];
				faceB.vertexColors[1] = hsl[3];
				faceB.vertexColors[2] = hsl[2];
			}

			surfaceGeometry.faces.push(faceA);
			surfaceGeometry.faces.push(faceB);		

			//create panel outlines.
			//We could outline each individual panel as a new line object, but that's super slow to render.
			//Alternatively, we could create just one line object and join up all the vertices - that sort
			//of works, and it's super fast, but it's ugly.
			//So a compromise is to have one line object for each chordwise strip of panels on the blade,
			//and one for each spanwise strip of wake panels.  Another way is to keep track of all the
			//vertices, and creatively use a three.LineSegments object.

			v = surfaceGeometry.vertices;		
			
			outlineGeom.vertices.push(v[k]);
			outlineGeom.vertices.push(v[k+1]);
			outlineGeom.vertices.push(v[k+1]);
			outlineGeom.vertices.push(v[k+2]);
			outlineGeom.vertices.push(v[k+2]);
			outlineGeom.vertices.push(v[k+3]);
			outlineGeom.vertices.push(v[k+3]);
			outlineGeom.vertices.push(v[k]);

			k+=4;
		}	
		
		surfaceGeometry.computeFaceNormals();

		var material = new THREE.MeshBasicMaterial( {
			color: surfaceColor,
			side: THREE.DoubleSide,
		});
		if(p.isScalarField){
			material.surfaceColor = 0xffffff;
			material.vertexColors = THREE.VertexColors;
		}	
		
		var surfaceMesh = new THREE.Mesh(surfaceGeometry, material);
		surfaceMesh.name = name_surfaceMesh;
		this.scene.add(surfaceMesh);
		
		surfaceGeometry.dispose();
		material.dispose();
		
		//mesh lines
		var lineMaterial = new THREE.LineBasicMaterial({
			color: 0x777777,
			linewidth: 0.1,
		});

		var meshLines = new THREE.LineSegments(outlineGeom, lineMaterial);
		meshLines.name = name_meshLines;
		this.scene.add(meshLines);
		
		outlineGeom.dispose();
		lineMaterial.dispose();
	}

	drawPanelCoordinateSystems(p){
		/*
		params = {
			mesh: model.s.mesh,
			isBladeMesh: true,
			isVisible = model.showBladePanelCoordinateSystems;
			name = "model.bladeCS";
		}
		*/

		//Create three Three.lineSegments objects, one for each panel axis.  Instead of having
		//one object per line (or even one object with three lines for each panel), we can draw
		//all lines that share the same material (like a1, a2 or n) as a single object.  Much faster.
		
		var name_a1 = p.name+".a1";
		var name_a2 = p.name+".a2";
		var name_n = p.name+".n";	
		
		this.removeFromScene(name_a1);
		this.removeFromScene(name_a2);	
		this.removeFromScene(name_n);
		
		if(!p.isVisible)return;
		
		var panelMesh;
		if(p.isBladeMesh)panelMesh = p.mesh.bladePanel;
		else panelMesh = p.mesh.wakePanel;
		
		var a1_geom = new THREE.Geometry();		
		var a2_geom = new THREE.Geometry();
		var n_geom = new THREE.Geometry();
		
		var p0, a1, a2, n;
		for(var i=0;i<panelMesh.length;i++){
			p0 = panelMesh[i].P[0].clone();
			a1 = panelMesh[i].a1.clone();
			a2 = panelMesh[i].a2.clone();			
			//scale n (unit normal) to be as long as a2 so it looks ok on the screen.			
			n = panelMesh[i].n.multiply(a2.NORM());
			
			a1 = a1.plus(p0);
			a2 = a2.plus(p0);
			n = n.plus(p0);

			//a1 direction
			a1_geom.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
			a1_geom.vertices.push(new THREE.Vector3(a1.x, a1.y, a1.z));		

			//a2 direction
			a2_geom.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
			a2_geom.vertices.push(new THREE.Vector3(a2.x, a2.y, a2.z));		

			//n direction
			n_geom.vertices.push(new THREE.Vector3(p0.x, p0.y, p0.z));
			n_geom.vertices.push(new THREE.Vector3(n.x, n.y, n.z));		
		}

		//create new Three.LineSegments objects.  Requires Threejs.72 or above.
		var material_a1 = new THREE.LineBasicMaterial({color: 0x0000ff});//blue for u direction
		var material_a2 = new THREE.LineBasicMaterial({color: 0x00ff00});//green for v direction
		var material_n = new THREE.LineBasicMaterial({color: 0xff0000});//red for surface normal		
		
		var a1_line = new THREE.LineSegments(a1_geom, material_a1);
		var a2_line = new THREE.LineSegments(a2_geom, material_a2);
		var n_line = new THREE.LineSegments(n_geom, material_n);
		
		a1_line.name = name_a1;
		a2_line.name = name_a2;
		n_line.name = name_n;
		
		this.scene.add(a1_line);
		this.scene.add(a2_line);
		this.scene.add(n_line);		
		
		a1_geom.dispose();
		a2_geom.dispose();
		n_geom.dispose();
		material_a1.dispose();
		material_a2.dispose();
		material_n.dispose();
	}

	drawCentreLine(p){
		/*
		params = {
			model:model.s, 
			isVisible: model.showInnerSurface, 
			name: "model.innerSurface.centreline", 
			surfaceIndex: i, 
			numUSteps:numUSteps, 
			numVSteps:numVSteps,
		};
		*/

		this.removeFromScene(p.name);
		if(!p.isVisible)return;

		//add a blue line to represent the zero-width surface at the centreline.
		var geom = new THREE.Geometry();
		
		var centre, centreCurves = p.model.curves[p.surfaceIndex];
		
		for(var i=0; i<centreCurves.length; i++){
			centre = centreCurves[i].c.P[0].p;
			geom.vertices.push(new THREE.Vector3(centre.x, centre.y, centre.z));
		}	
		
		var material = new THREE.LineBasicMaterial({color: 0x0000ff});

		var line = new THREE.Line(geom, material);
		line.name = p.name;
		
		this.scene.add(line);

		material.dispose();
		geom.dispose();
	}

	drawAxes(p){
		/*
		params = {
			model: model,
			name: "model.axes",
			isVisible: model.showAxes,
		}
		*/

		var name_x = p.name+".x";
		var name_y = p.name+".y";	
		var name_z = p.name+".z";
		
		var name_plane_xy = p.name+".xy";
		var name_plane_xz = p.name+".xz";
		var name_plane_yz = p.name+".yz";

		this.removeFromScene(name_x);
		this.removeFromScene(name_y);	
		this.removeFromScene(name_z);
		this.removeFromScene(name_plane_xy);
		this.removeFromScene(name_plane_xz);
		this.removeFromScene(name_plane_yz);
		
		if(!p.isVisible)return;

		var surface = p.model;	
		
		//create and add global origin axes
		var length = surface.s.curves[0][surface.s.curves[0].length-1].c.P[0].p.minus(surface.s.curves[0][0].c.P[0].p).NORM()/10;

		var geom, material;
		//x-axis	
		material = new THREE.LineBasicMaterial({color: 0xff0000});	
		geom = new THREE.Geometry();		
		
		geom.vertices.push(new THREE.Vector3(0., 0., 0.));
		geom.vertices.push(new THREE.Vector3(length, 0., 0.));		
		
		var line_x = new THREE.Line(geom, material);
		line_x.name = name_x;
		this.scene.add(line_x);

		geom.dispose();
		material.dispose();

		//y-axis
		material = new THREE.LineBasicMaterial({color: 0x00ff00});
		geom = new THREE.Geometry();		
		
		geom.vertices.push(new THREE.Vector3(0., 0., 0.));
		geom.vertices.push(new THREE.Vector3(0., length, 0.));		
		
		var line_y = new THREE.Line(geom, material);
		line_y.name = name_y;
		this.scene.add(line_y);

		geom.dispose();
		material.dispose();
		
		//z-axis
		material = new THREE.LineBasicMaterial({color: 0x0000ff});
		geom = new THREE.Geometry();		
		
		geom.vertices.push(new THREE.Vector3(0., 0., 0.));
		geom.vertices.push(new THREE.Vector3(0., 0., length));		
		
		var line_z = new THREE.Line(geom, material);
		line_z.name = name_z;
		this.scene.add(line_z);

		geom.dispose();
		material.dispose();

		//Planes.	
		if(!p.hideAxisPlanes){
			//find max and min coordinates of the 3D object.
			var x0=1.E7;
			var x1=-1.E7;
			var y0=1.E7;
			var y1=-1.E7;
			var z0=1.E7;
			var z1=-1.E7;
			
			var j;
			for(var i=0;i<surface.s.curves[0].length;i++){
				for(j=0;j<surface.s.curves[0][i].c.P.length;j++){
					if(surface.s.curves[0][i].c.P[j].p.x<x0)x0 = surface.s.curves[0][i].c.P[j].p.x;
					else if(surface.s.curves[0][i].c.P[j].p.x>x1)x1 = surface.s.curves[0][i].c.P[j].p.x;
					
					if(surface.s.curves[0][i].c.P[j].p.y<y0)y0 = surface.s.curves[0][i].c.P[j].p.y;
					else if(surface.s.curves[0][i].c.P[j].p.y>y1)y1 = surface.s.curves[0][i].c.P[j].p.y;
					
					if(surface.s.curves[0][i].c.P[j].p.z<z0)z0 = surface.s.curves[0][i].c.P[j].p.z;
					else if(surface.s.curves[0][i].c.P[j].p.z>z1)z1 = surface.s.curves[0][i].c.P[j].p.z;
				}
			}
			
			var margin_x = Math.abs(x1-x0)*0.3;
			var margin_y = Math.abs(y1-y0)*0.3;		
			var margin_z = Math.abs(z1-z0)*0.3;
			var margin = Math.min(margin_x, margin_y, margin_z);
			var corners = {
				x: [x0-margin,x1+margin],
				y: [y0-margin,y1+margin],
				z: [z0-margin,z1+margin]
			};
			
			//////
			//xy
			this.drawAxisPlane({x:1, y:1, z:0, c:corners, name:name_plane_xy, isVisible:p.isVisible});
			//////
			//xz
			this.drawAxisPlane({x:1, y:0, z:1, c:corners, name:name_plane_xz, isVisible:p.isVisible});
			//////
			//yz
			this.drawAxisPlane({x:0, y:1, z:1, c:corners, name:name_plane_yz, isVisible:p.isVisible});
		}	
	}

	drawAxisPlane(p){
		/*
		params = {x:1, y:1, z:0, c:corners, name:p.name+".xy", isVisible:p.isVisible}
		*/	

		this.removeFromScene(name);
		if(!p.isVisible)return;

		//create and add global origin axes	
		var geom;
		
		var name_line = p.name+".line";
		var name_plane = p.name+".plane";
		
		var geom = new THREE.Geometry();
		
		var planeSideWidth, planeSideHeight;
		
		if(p.x==1 && p.y==1){
			geom.vertices.push(new THREE.Vector3(p.c.x[0], p.c.y[0], 0));
			geom.vertices.push(new THREE.Vector3(p.c.x[1], p.c.y[0], 0));		
			geom.vertices.push(new THREE.Vector3(p.c.x[1], p.c.y[1], 0));
			geom.vertices.push(new THREE.Vector3(p.c.x[0], p.c.y[1], 0));	
			geom.vertices.push(new THREE.Vector3(p.c.x[0], p.c.y[0], 0));
			planeSideWidth = p.c.x[1] - p.c.x[0];
			planeSideHeight = p.c.y[1] - p.c.y[0];
		}	
		else if(p.x==1 && p.z==1){
			geom.vertices.push(new THREE.Vector3(p.c.x[0], 0, p.c.z[0]));
			geom.vertices.push(new THREE.Vector3(p.c.x[1], 0, p.c.z[0]));		
			geom.vertices.push(new THREE.Vector3(p.c.x[1], 0, p.c.z[1]));
			geom.vertices.push(new THREE.Vector3(p.c.x[0], 0, p.c.z[1]));	
			geom.vertices.push(new THREE.Vector3(p.c.x[0], 0, p.c.z[0]));
			planeSideHeight = p.c.x[1] - p.c.x[0];
			planeSideWidth = p.c.z[1] - p.c.z[0];			
		}	
		else if(p.y==1 && p.z==1){
			geom.vertices.push(new THREE.Vector3(0, p.c.y[0], p.c.z[0]));
			geom.vertices.push(new THREE.Vector3(0, p.c.y[1], p.c.z[0]));		
			geom.vertices.push(new THREE.Vector3(0, p.c.y[1], p.c.z[1]));
			geom.vertices.push(new THREE.Vector3(0, p.c.y[0], p.c.z[1]));	
			geom.vertices.push(new THREE.Vector3(0, p.c.y[0], p.c.z[0]));		
			planeSideWidth = p.c.y[1] - p.c.y[0];
			planeSideHeight = p.c.z[1] - p.c.z[0];			
		}	
		
		var materialOutline = new THREE.LineBasicMaterial({color: 0xffff00});
		
		var container = new THREE.Object3D();
		container.name = p.name;
		
		//plane outline
		var line = new THREE.Line(geom, materialOutline);
		line.name = name_line;
		line.visible = p.isVisible;	
		container.add(line);

		geom.dispose();
		materialOutline.dispose();
		
		//////Shaded transparent planes that fit inside the outlines.
		var planeGeometry = new THREE.PlaneBufferGeometry( planeSideWidth, planeSideHeight,1,1);
		
		var material = new THREE.MeshBasicMaterial( {
			color: 0xffff00, //yellow
			side: THREE.DoubleSide,
			transparent: true, 
			opacity: 0.04,
		});				

		var plane = new THREE.Mesh(planeGeometry, material);
		plane.name = name_plane;
		plane.visible = p.isVisible;	
		container.add(plane);
		
		planeGeometry.dispose();
		material.dispose();

		if(p.x==1 && p.y==1){
			plane.rotation.x = 0.;
			plane.rotation.y = 0.;
			plane.rotation.z = 0.;
			
			plane.position.x = p.c.x[0] + (p.c.x[1]-p.c.x[0])/2;
			plane.position.y = p.c.y[0] + (p.c.y[1]-p.c.y[0])/2;
			plane.position.z = 0;
		}	
		else if(p.x==1 && p.z==1){
			plane.rotation.x = Math.PI/2.;
			plane.rotation.y = 0.;
			plane.rotation.z = Math.PI/2.;

			
			plane.position.x = p.c.x[0] + (p.c.x[1]-p.c.x[0])/2;
			plane.position.y = 0.;
			plane.position.z = p.c.z[0] + (p.c.z[1]-p.c.z[0])/2;
		}	
		else if(p.y==1 && p.z==1){
			plane.rotation.x = 0.;
			plane.rotation.y = Math.PI/2.;
			plane.rotation.z = Math.PI/2.;

			plane.position.x = 0;
			plane.position.y = p.c.y[0] + (p.c.y[1]-p.c.y[0])/2;
			plane.position.z = p.c.z[0] + (p.c.z[1]-p.c.z[0])/2;
		}
		
		this.scene.add(container);
	}
	//end of DrawModel methods
}//end of DrawModel
//
//
Number.prototype.round = function(places) {
  var result = +(Math.round(this + "e+" + places)  + "e-" + places);
  if(isNaN(result))return 0.;
  return result;
}	