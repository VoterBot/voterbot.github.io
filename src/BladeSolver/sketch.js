"use strict";
/******
Copyright 2000-2016 Mark Hampsey // escapee.org

Portions of this code are adapted from routines in 'The NURBS Book' by Piegl and Tiller:
http://www.springer.com/gp/book/9783642973857
*******/
function loadSketch(){
	var canvas_width = 800;
	var canvas_height = 500;

	////////
	//The controlling object. Contains current SolidObject and screen elements.
	var solidModel = new SolidObject();	//sets the default object: a B-spline model of a 5 kW blade.
	
	//set up tab with blade/wake mesh.
	var canvasDiv = document.getElementById("drawingCanvas");
	solidModel.draw.bladeAndWake = new DrawModel(canvas_width, canvas_height, canvasDiv);

	//set up blade design tab.
	var bladeDesignCanvasDiv = document.getElementById("bladeDesignCanvas");
	solidModel.draw.bladeDesign = new DrawModel(canvas_width, 200, bladeDesignCanvasDiv, true);

	//set up the text DOM element for area and volume output.
	solidModel.outputToScreen.areaAndVolume = document.getElementById("areaAndVolumeOutput");
	solidModel.outputToScreen.areaAndVolume.style.left = "20px";
	solidModel.outputToScreen.areaAndVolume.style.width = canvas_width+"px";		
	//solidModel.outputToScreen.areaAndVolume.style.top = (canvas_height-35)+"px";
	solidModel.outputToScreen.areaAndVolume.style.top = "40px";
	
	//set up the text DOM element for solver results output.
	var solverDiv = document.getElementById("solverResultsOutput");
	solverDiv.style.display = "none";		
	solverDiv.style.left = (canvas_width-245)+"px";
	solverDiv.style.top = (canvas_height-290)+"px";
	solverDiv.style.width = "230px";		
	solverDiv.style.height = "290px";
	solverDiv.style.backgroundColor = "lightblue";		
	solverDiv.style.opacity = "0.75";
	solverDiv.style.color = "rgba(0, 0, 0, 1.0)";
	solverDiv.style.fontSize = "12px";
	solverDiv.style.fontFamily = "sans-serif";		
	solverDiv.style.padding = "5px";
	
	solidModel.outputToScreen.solverResults = solverDiv;	

	//
	var deleteButton = document.getElementById('solverResultsDeleteButton');
	deleteButton.style.top = solverDiv.style.top;
	deleteButton.style.left = (parseInt(solverDiv.style.left)+parseInt(solverDiv.style.width)-45)+"px";
	
	solidModel.outputToScreen.solverResultsButtons.deleteButton = deleteButton;
	solidModel.outputToScreen.solverResultsButtons.deleteButton.onclick = function(){
		solidModel.solutions.remove();
	};		
	
	//
	var saveButton = document.getElementById('solverResultsSaveButton');
	saveButton.style.top = solverDiv.style.top;
	saveButton.style.left = (parseInt(deleteButton.style.left)+18)+"px";		
	
	solidModel.outputToScreen.solverResultsButtons.saveButton = saveButton;
	solidModel.outputToScreen.solverResultsButtons.saveButton.onclick = function(){
		solidModel.solutions.compressAndSave();
	};		
	
	//
	var closeButton = document.getElementById('solverResultsCloseButton');
	closeButton.style.top = solverDiv.style.top;
	closeButton.style.left = (parseInt(saveButton.style.left)+18)+"px";
	
	solidModel.outputToScreen.solverResultsButtons.closeButton = closeButton;
	solidModel.outputToScreen.solverResultsButtons.closeButton.onclick = function(){
		solidModel.outputToScreen.showSolverResults(false);
	};		

	//draw everything.
	solidModel.refresh();
	
	window.onbeforeunload = function(){
		solidModel.destroy();
	}	
};
////
//SolidObject class.  Contains a NURBS model, as well as display meshes and panel solver meshes for it.
class SolidObject{
	constructor(model){
		if(!model)model = this.getDefaultModel();
		
		this.s = model.NURBS_model; //the NURBS model.
		
		//display flags
		this.showAllCurves = false;
		this.showAxes = false;
		this.showOuterSurface = false;
		this.showInnerSurface = false;
		this.showSurfaceDerivatives = false;
		//
		this.showBladePanelMesh = true;
		this.showWakePanelMesh = true;
		this.showBladePanelCoordinateSystems = false;
		this.showWakePanelCoordinateSystems = false;
		//
		this.showResults;
		
		//"Model" and "Mesh" menus.
		if(model.modelSettings){
			this.modelSettings = model.modelSettings;
			this.modelSettings.settings.theObject = this.s;
			this.modelSettings.settings.parent = this;
		}	
		if(model.meshSettings){
			this.meshSettings = model.meshSettings;
			this.meshSettings.settings.theObject = this.s.mesh;
			this.meshSettings.settings.parent = this;
		}
		
		//"Display" menu.  Doesn't need to be a "Settings" object, because it never needs to be reset.
		var self = this;
		this.displaySettings = {
			//model
			showOuterSurface: {value: self.showOuterSurface, name: "outer surface"},			
			showInnerSurface: {value: self.showInnerSurface, name: "inner isosurfaces"},
			showSurfaceDerivatives: {value: self.showSurfaceDerivatives, name: "derivatives"},
			showAllCurves: {value: self.showAllCurves, name: "all curves"},
			showAxes: {value: self.showAxes, name: "show axes"},
			//
			folderName: "Display",			
			theObject: self,
		};
		if(model.meshSettings){
			//mesh
			this.displaySettings.showBladePanelMesh = {value: this.showBladePanelMesh, name: "Blade panels"};
			this.displaySettings.showWakePanelMesh = {value: this.showWakePanelMesh, name: "Wake panels"};
			this.displaySettings.showBladePanelCoordinateSystems = {value: this.showBladePanelCoordinateSystems, name: "Blade CS"};
			this.displaySettings.showWakePanelCoordinateSystems = {value: this.showWakePanelCoordinateSystems, name: "Wake CS"};
		}
		else{
			//disable the mesh menu options.
			this.showBladePanelMesh;
			this.showWakePanelMesh;
			this.showBladePanelCoordinateSystems;
			this.showWakePanelCoordinateSystems;
			//enable curve display.
			this.showAllCurves = true;
		}
		//
		
		//"Results" menu.  Doesn't need to be a "Settings" object, because it never gets reset.
		this.resultsSettings = {
			//results.
			set_Cp: {name: "Cp"},		
			set_phi: {name: "phi"},
			set_dphi_dn: {name: "dphi_dn"},
			set_pressure: {name: "pressure"},
			set_x_on_C: {name: "x_on_C"},
			set_r_on_R: {name: "r_on_R"},
			//
			folderName: "Results",		
			theObject: self,
		};
		this.set_Cp = function(){this.showResults = {field:"Cp", reverseYAxis:true};};
		this.set_phi = function(){this.showResults = {field:"phi", units:"m"+String.fromCharCode(178)+"/s"};};
		this.set_dphi_dn = function(){this.showResults = {field:"dphi_dn", units:"m/s"};};
		this.set_pressure = function(){this.showResults = {field:"pressure", units:"Pascals", reverseYAxis:true};};	
		this.set_x_on_C = function(){this.showResults = {field:"x_on_C", hidePlots:true};};
		this.set_r_on_R = function(){this.showResults = {field:"r_on_R", hidePlots:true};};	
		//
		
		/////
		//
		this.draw = {}; //a container for THREE.js canvas drawing objects.
		//
		this.areaAndVolumeWorker; //a web worker that calculates things in the background.
		//
		this.solverWorker; //the web worker for the panel solver.
		//
		this.solutions = new Solutions(this); //a collection of solved panel meshes.
		////
		//the dat.gui menu.
		this.solverMenu = {gui: undefined, controlDiv:undefined};
		this.designMenu = {gui: undefined, controlDiv:undefined};	
		this.resultsMenu = {gui: undefined, controlDiv:undefined};

		//an object for text overlays which show progress of solvers etc.
		this.outputToScreen = {
			areaAndVolume: undefined,
			//
			solverResults: undefined,
			solverResultsButtons: {
				deleteButton: undefined,
				saveButton: undefined,
				closeButton: undefined,
			},	
			showSolverResults: function(isShow, hasButtons){
				if(isShow)this.solverResults.style.display = "inline";
				else this.solverResults.style.display = "none";
			
				for(var key in this.solverResultsButtons){
					if(hasButtons)this.solverResultsButtons[key].style.display = "inline-block";			
					else this.solverResultsButtons[key].style.display = "none";
				}	
			},
		}
		////
	}	
	
	//methods
	resetModel(skipReset){
		//reset NURBS model back to default settings.
		if(this.modelSettings && !skipReset)this.s.setBladeSettings(this.modelSettings.getValues(true));
		//redraw everything.
		this.refresh({model:true, chordAndTwist:{isReset:true}, gui:true});
	}

	resetPanelMesh(skipReset){
		//retrieve default panel mesh settings.
		var meshSettings = undefined;
		if(this.meshSettings && !skipReset)meshSettings = this.meshSettings.getValues(true);
		//redraw mesh.
		this.refresh({mesh:true, meshDefaults:meshSettings, gui:true});	
	}

	plotChordAndTwist(){
		//create interactive chord and twist plots.		
		var parameters = {
			plotContainer: document.getElementById('designPlotContainer'),
			plotWidth: 350,
			plotMargin: {left:50, right:35, top:35, bottom:50},
			numPlotsAcrossPage: 2,
			hasBorder: false,
		};	

		//remove any existing plots.
		while(parameters.plotContainer.hasChildNodes()){  
			parameters.plotContainer.removeChild(parameters.plotContainer.firstChild);
		}			
		
		var plotPanel = new PlotPanel(parameters);
		
		//
		//// add plots
		//

		var bladeDesign = this.s.getChordAndTwistDistribution(); //the desired chord and twist distributions.
		var current = this.s.getCurrentChordAndTwist(); //the actual chord and twist of the current blade.

		var chord = [], twist = [], p;
		for(var i=0;i<bladeDesign.template.length;i++){
			p = bladeDesign.template[i];
			chord.push(new WeightedPoint_2D(p.p.x, p.p.y, p.w));
			twist.push(new WeightedPoint_2D(p.p.x, p.p.z, p.w));				
		}
		
		var self = this;
		var chordPlot = {
			series: [chord, current.chordDistribution],
			bsplineCurveFit: bladeDesign.bspline.extract2DCurve({isXY:true}),
			//xTicStep:0.1, 
			title:"chord", 
			xAxisLabel:"r/R", 
			yAxisLabel:"normalised chord (chord/R)",
			//hasBorder: true,
			callback: function(data){
				var curve = bladeDesign.bspline; //directly modify the curve/twist definition object of the model.

				for(var i=0;i<this.bsplineCurveFit.P.length;i++){
					curve.c.P[i].p.x = this.bsplineCurveFit.P[i].p.x;
					curve.c.P[i].p.y = this.bsplineCurveFit.P[i].p.y;						
				}

				//console.log("UPDATING DISPLAY LIKE THIS CAUSES BIG MEMORY LEAKS.  FIX.");					
				self.refresh({model:true});
			}
		};
		plotPanel.addPlot(chordPlot);

		var twistPlot = {
			series: [twist, current.twistDistribution],
			bsplineCurveFit: bladeDesign.bspline.extract2DCurve({isXZ:true}),
			bsplineAllowDrag: {x:false, y:true},
			//xTicStep:0.1, 
			title:"twist", 
			xAxisLabel:"r/R", 
			yAxisLabel:"twist (deg)",
			//hasBorder: true,
			callback: function(data){
				var curve = bladeDesign.bspline;

				for(var i=0;i<this.bsplineCurveFit.P.length;i++){
					curve.c.P[i].p.x = this.bsplineCurveFit.P[i].p.x;
					curve.c.P[i].p.z = this.bsplineCurveFit.P[i].p.y;
				}
				
				self.refresh({model:true});
			}
		};
		plotPanel.addPlot(twistPlot);			
	}

	plotResults(){
		//create plots of spanwise scalar fields.
		var parameters = {
			plotContainer: document.getElementById('resultsPlotContainer'),
			plotWidth: 350,
			plotMargin: {left:50, right:35, top:35, bottom:50},
			numPlotsAcrossPage: 2,
			hasBorder: false,
		};	

		//remove any existing plots.
		while(parameters.plotContainer.hasChildNodes()){  
			parameters.plotContainer.removeChild(parameters.plotContainer.firstChild);
		}			
		
		if(this.showResults.hidePlots)return;
		
		var plotPanel = new PlotPanel(parameters);
		
		//
		//// add plots
		//

		var self = this;
		
		var i, j, fieldPlot, series;
		var field = this.showResults.field, units = this.showResults.units, reverseYAxis = this.showResults.reverseYAxis;
		
		if(units === undefined)units = "";
		else units = " ("+units+")";
		
		var panels = this.s.mesh.bladePanel;	

		var index = 0;
		for(var j=0;j<this.s.mesh.mPanels;j++){
			series = [];
			for(i=0;i<this.s.mesh.nPanels;i++){
				series.push(new Point_2D(panels[index].x_on_C, panels[index][field]));
				index++;
			}
			
			fieldPlot = {
				series: [series],
				title:"r/R = "+panels[index-1].r_on_R.round(4), 
				xAxisLabel:"x/C", 
				yAxisLabel:field+units,
				reverseYAxis: reverseYAxis, //suction side up for pressure plots.
			};
			plotPanel.addPlot(fieldPlot);
		}	
		
		//console.log(this.s.mesh.bladePanel);	
	}

	outputAreaAndVolume(){
		this.outputToScreen.areaAndVolume.innerHTML = "Surface area: calculating ...";
		this.outputToScreen.areaAndVolume.innerHTML += "&nbsp;&nbsp;&nbsp;Volume: calculating ...";

		///////
		//If we've interrupted a currently-running worker, kill it.
		if(this.areaAndVolumeWorker){
			this.areaAndVolumeWorker.terminate();
			this.areaAndVolumeWorker = undefined;
		}

		//Make a new web worker.		
		var worker = new Worker('NURBSSolid.js');
		this.areaAndVolumeWorker = worker;
		
		var area = undefined, volume = undefined;
		var output = this.outputToScreen.areaAndVolume;
		output.set = function(text, isClear){
			if(isClear)this.innerHTML = "";
			this.innerHTML += text;
		};
		
		worker.onmessage = function(e){
			if(e.data.result.area)area = e.data.result.area.round(5);
			if(e.data.result.volume)volume = e.data.result.volume.round(5);
			
			//console.log(area+" : "+volume);

			if(area)output.set("Surface area: "+area+" m<sup>2</sup>", true);
			else output.set("Surface area: calculating ...", true);
			if(volume)output.set("&nbsp;&nbsp;&nbsp;Volume: "+volume+" m<sup>3</sup>", false);
			else output.set("&nbsp;&nbsp;&nbsp;Volume: calculating ...", false);

			if(e.data.message){
				//console.log(e.data.message);
			}

			if(area && volume){
				worker.terminate();
				worker = undefined;
			}
		};
		
		worker.onerror = function(e) {
			console.log('Error: Line ' + e.lineno + ' in ' + e.filename + ': ' + e.message);
		};
		
		//start the worker
		worker.postMessage({'cmd': 'areaAndVolume', 'solidModel': this.s,});
	}

	stopSolver(){
		this.solverWorker = this.s.mesh.stopSolver(this.solverWorker);		
	}

	solve(){
		//remove any existing results tabs.
		this.setupResultsTab({isHide:true});	

		//toggle the solver on and off.
		if(this.solverWorker){
			this.outputToScreen.showSolverResults(false);
			this.stopSolver();
		}	
		else{
			console.log("STARTING SOLVER ...");
			this.outputToScreen.showSolverResults(true, false);

			//an object with functions to handle messages from the web worker.
			var self = this;
			var resultsHandler = {
				updateMeshDisplay: function(x){
					self.outputToScreen.showSolverResults(true, true); //enable results buttons.
					self.s.mesh = x; //update the panel mesh.
					self.solutions.add(self.s, self.outputToScreen.solverResults.innerHTML);
					self.stopSolver(); //remove the solver web worker.
					
					//redraw.
					self.refresh({gui:true, meshRedraw:true, results:true});
				},
				
				setSolverResults: function(text){
					self.outputToScreen.solverResults.innerHTML = text;
				},	
			};
			
			//launch the solver web worker.
			this.solverWorker = this.s.mesh.solve(resultsHandler);
		}
	}

	setupGUI(){
		var self = this; //the solidModel object.
		
		//remove any existing menus.
		if(self.solverMenu.gui){
			//remove the old solver gui
			self.solverMenu.controlDiv.removeChild(self.solverMenu.gui);
			self.solverMenu = {gui: undefined, controlDiv:undefined};
		}	
		//
		if(self.designMenu.gui){
			//remove the old design gui
			self.designMenu.controlDiv.removeChild(self.designMenu.gui);
			self.designMenu = {gui: undefined, controlDiv:undefined};
		}
		//
		if(self.resultsMenu.gui){
			//remove the old results gui	
			self.resultsMenu.controlDiv.removeChild(self.resultsMenu.gui);
			self.resultsMenu = {gui: undefined, controlDiv:undefined};
		}	

		//ensure pointers are referencing updated objects.
		if(this.displaySettings)this.displaySettings.theObject = this;
		if(this.resultsSettings)this.resultsSettings.theObject = this;	
		if(this.modelSettings){
			this.modelSettings.settings.theObject = this.s;
			this.modelSettings.settings.parent = this;
		}
		if(this.meshSettings){
			this.meshSettings.settings.theObject = this.s.mesh;
			this.meshSettings.settings.parent = this;
		}	
		
		//create menus.
		var gui_solver = new dat.GUI({autoPlace: false});
		var gui_design = new dat.GUI({autoPlace: false});	
		var gui_results = new dat.GUI({autoPlace: false});
		
		//a function to create a menu from one of the parameter lists of a SolidObject.
		function setModelMenu(settings, gui){
			var modelFolder = gui.addFolder(settings.folderName);
			modelFolder.close();
			//
			var key, menuItem, theObject, isReset;
			for(key in settings){
				isReset = false;
				if(key.indexOf("reset")==0)isReset = true;
				
				if(isReset)theObject = settings.parent;
				else theObject = settings.theObject;
				
				if((theObject[key] !== undefined && settings[key].name && !settings[key].isHidden) || isReset){
					menuItem = modelFolder.add(theObject, key).name(settings[key].name);
					if(settings[key].isInteger)menuItem.step(1);
					if(settings[key].min != undefined)menuItem.min(settings[key].min);

					menuItem.onFinishChange(function(value){
						if(!this.property.startsWith("reset")){
							//redraw everything.  "reset" has its own refresh routine.
							var refresh = {gui:false, model:false, mesh:false, results:false, display:false};
							 
							if(theObject.curves)refresh.model = true; //are we changing the Bspline object?
							else if(theObject.wakePanel)refresh.mesh = true;
							else if(settings.folderName == "Display")refresh.display = this.property;
							else if(settings.folderName == "Results")refresh = {results:true, gui:true};
							
							self.refresh(refresh);
						}	
					});	
				}	
			}
		}

		//create one menu folder each for: NURBSSolid, its panel mesh and the display choices.
		if(self.modelSettings)setModelMenu(self.modelSettings.settings, gui_design); //NURBSSolid.
		if(self.meshSettings)setModelMenu(self.meshSettings.settings, gui_solver); //Panel mesh.
		setModelMenu(self.displaySettings, gui_solver); //Display choices.	
		if(self.draw.results)setModelMenu(self.resultsSettings, gui_results); //Panel mesh.	

		//add a solve button		
		if(self.meshSettings){	
			var solveButton = gui_solver.add(self, "solve");
			solveButton.onFinishChange(function(value){
				for(var x in gui_solver.__folders)gui_solver.__folders[x].close(); //hide all open folders.
				if(self.solverWorker)this.name("STOP SOLVER");
				else this.name("solve");
			});

			//add a folder for solutions.
			var solutionsFolder = gui_solver.addFolder("Solutions");
			solutionsFolder.close();
			
			//add a list of solutions, if they exist.
			for(var key in self.solutions.model){
				solutionsFolder.add(self.solutions.model[key], "set").name(self.solutions.model[key].name);
			}
			
			//add a 'load' button.
			var loadSolutionButton = solutionsFolder.add(self.solutions, "loadFromFile").name("load saved");
		}	
		
		//retrieve the DOM container for the solver menu, and add the gui object to it.
		self.solverMenu.controlDiv = document.getElementById("guiControls_solver");  
		self.solverMenu.controlDiv.style.left = (this.draw.bladeAndWake.canvas_width-gui_solver.width)+"px";
		self.solverMenu.gui = self.solverMenu.controlDiv.appendChild(gui_solver.domElement);
		//
		//retrieve the DOM container for the design menu, and add the gui object to it.
		self.designMenu.controlDiv = document.getElementById("guiControls_design");  
		self.designMenu.controlDiv.style.left = (this.draw.bladeAndWake.canvas_width-gui_design.width)+"px";
		self.designMenu.gui = self.designMenu.controlDiv.appendChild(gui_design.domElement);
		//
		//retrieve the DOM container for the results menu, and add the gui object to it.
		if(this.draw.results){
			self.resultsMenu.controlDiv = document.getElementById("guiControls_results");  
			self.resultsMenu.controlDiv.style.left = (this.draw.results.canvas_width-gui_results.width-50)+"px";
			self.resultsMenu.gui = self.resultsMenu.controlDiv.appendChild(gui_results.domElement);
		}	
	}//end setupGUI

	setupResultsTab(flags){
		if(flags && flags.isHide && this.draw.results === undefined)return; //already hidden.

		//hide/unhide the 'results' tab.
		if(flags && flags.isHide){
			 //hide the tab.
			document.getElementById("resultsTabLabel").style.display = "none";
			return;
		}	
		else{
			//reveal the tab.
			document.getElementById("resultsTabLabel").style.display = "inline-block";
			
			//initialise the results canvas object if it doesn't already exist.
			if(!this.draw.results){
				//retrieve DOM container for 3D canvas.
				var resultsModelDiv = document.getElementById("resultsModelContainer");
			
				this.draw.results = new DrawModel(this.draw.bladeDesign.canvas_width, this.draw.bladeDesign.canvas_height, resultsModelDiv, true);
			}	
			
			//retrieve and hide the gradient canvas - container for a nice colour bar for scalar field legends.
			var gradientCanvas = document.getElementById("resultsGradientCanvas");
			gradientCanvas.style.display = "none"; //hide it, initially.
			
			//draw the 3D blade mesh, with a colour gradient overlay for various solver results.
			this.draw.results.addSolverResultsToScene(this, gradientCanvas);	
			
			//console.log(this.solutions.model[this.solutions.currentSolutionKey]);
			
			//create results plots.
			this.plotResults();
		}	
	}

	getDefaultModel(){
		//load the two bspline aerofoils from iframes embedded in the page.
		var baseAerofoil = new SurfaceCurve();
		baseAerofoil.initialiseFromFile(this.loadData("7062_bspline"));
		var hubAerofoil = new SurfaceCurve();
		hubAerofoil.initialiseFromFile(this.loadData("rectangle_bspline"));
		
		///////
		//NURBSSolid
		//create the b-spline blade (an extension of a NURBSSolid).
		var bladeSettings = new Settings(
			{
				numFoils: {value: 15, name:"Number aerofoils", min:4, isInteger: true},
				rBlade: {value: 2.5, name: "Blade length (m)", min:0.1},
				pitchAngle: {value: 0., name: "Pitch angle (deg)"},
				transitionPercent: {value: 28.8, name: "Transition \%", min:10.},
				modelTransition: {value: true, name: "Model transition"},
				//
				baseFoil: {value: baseAerofoil},
				hubFoil: {value: hubAerofoil},
				//
				folderName: "Blade parameters",			
				resetModel: {name:"reset"},
			}
		);	

		var bsplineBlade = new BsplineBlade();
		bsplineBlade.createBlade(bladeSettings.getValues());
		
		///////
		//Panel Mesh
		//create the panel method surface mesh based on the b-spline blade.
		var meshSettings = new Settings(
			{
				numBlades: {value: 1, name:"Number of blades", min:1, isInteger:true},
				nPanels: {value: 10, name:"nPanels", min:4, isInteger:true},
				mPanels: {value: 10, name:"mPanels", min:2, isInteger:true},
				wPanels: {value: 10, name:"wPanels", min:4, isInteger:true},			
				turns: {value: 10, name:"wake turns", min:1, isInteger:true},
				lambda: {value: 8.5, name: "Tip Speed Ratio"},
				U_0: {value: 10, name: "Wind speed (m/s)"},			
				isHelicalWake: {value: true, name: "Helical wake?", isHidden:true},
				allowWakeRollup: {value: false, name: "Wake rollup?"},			
				resetPanelMesh: {name:"reset"},
				//
				folderName: "Panel mesh parameters",			
			}
		);	
		bsplineBlade.mesh = new BladeMesh(bsplineBlade, meshSettings.getValues());
		
		//default mesh parameters.
		return {NURBS_model: bsplineBlade, modelSettings: bladeSettings, meshSettings: meshSettings};
	}

	refresh(flags){
		//refresh menus, plots and anything else that changes when the model changes.
		//flags: {gui:bool, model:bool, modelRedraw:bool, mesh:bool, chordAndTwist:LeastSquaresCurve()}

		//if model settings have changed, recreate it.
		if(flags && flags.model)
			this.s.initialiseModel(flags.chordAndTwist); //recreate the solid object, with new design if available.
			
		//create a new default mesh if the geometry has changed.
		if(flags && (flags.model || flags.mesh) && this.s.mesh)
			this.s.mesh.setup(this.s, flags.meshDefaults);	

		//if model has changed, or on first-run, output things concerning new model.
		if(!flags || flags.model || flags.modelRedraw){
			//plot chord and twist of blade.
			this.plotChordAndTwist();
			
			//calculate and output surface area and volume for this model.
			this.outputAreaAndVolume();
		}

		var displayHasChanged = undefined;
		if(flags && flags.display)displayHasChanged = flags.display;

		//redraw blade/mesh/wake on solver tab.
		if(!flags || flags.model || flags.modelRedraw || flags.mesh || flags.meshRedraw || displayHasChanged)
			this.draw.bladeAndWake.addModelToScene(this, displayHasChanged); 
		
		//redraw blade design tab.
		if(!flags || flags.model || flags.modelRedraw)
			this.draw.bladeDesign.addBladeDesignToScene(this);	
		
		//reveal/hide the results tab - model changes invalidate previous solution.
		if(flags && flags.results)this.setupResultsTab();
		else this.setupResultsTab({isHide:true});
		
		//Refresh/collapse the menus if anything has changed.
		if(!flags || flags.gui)this.setupGUI();
	}

	loadData(fileID){
			//extract data from the datafile loaded into an iframe in the HTML.
			//NOTE:  To run locally, Chrome requires this switch on the command line: 
			//  --allow-file-access-from-files
			//Otherwise, it gets angry at you trying to load data from a local iframe.
			//Should work fine when served by a web server.
			var oFrame = document.getElementById(fileID);
			return oFrame.contentWindow.document.body.childNodes[0].innerHTML;
	}

	destroy(){
		console.log("DESTROY!!");
	}
	//end of SolidObject methods
}//end of SolidObject
////
////
//an object to store solutions.
class Solutions{
	constructor(parent){
		this.parent = parent;
		//a collection of BsplineBlades and their solved meshes.
		this.model = {}; //container for new solutions.
		this.currentSolutionKey = undefined;	
	}
	
	//methods
	add(x, resultsText){
		//add a solved model to the collection.	
		var self = this;
		var key = Date.now(); //tag each solution with the current time in milliseconds.
		self.model[key] = {
			model: x.clone(),
			resultsText: resultsText,
			name: x.mesh.numBlades+"b|"+x.mesh.nPanels+"n|"+x.mesh.mPanels+"m|"+x.mesh.wPanels+"w",
			set: function(){
				self.set(key);
			},	
		};
		this.currentSolutionKey = key;
	}

	set(key){
		//set the model in the current display object.
		this.parent.s = null;
		this.parent.s = this.model[key].model.clone();
		this.currentSolutionKey = key;
		this.updateDisplay(key);
	}

	remove(){
		//remove the currently-selected solution.
		if(this.currentSolutionKey){
			delete this.model[this.currentSolutionKey];
			this.currentSolutionKey = undefined;
			this.updateDisplay();
		}	
	}
	
	saveToFile(textToWrite){
		//example from: 
		//https://thiscouldbebetter.wordpress.com/2012/12/18/loading-editing-and-saving-a-text-file-in-html5-using-javascrip/
		var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
		var fileNameToSaveAs = this.currentSolutionKey+".dat";			
		var downloadLink = document.createElement("a");
		downloadLink.download = fileNameToSaveAs;
		downloadLink.innerHTML = "Download File";
		
		downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
		downloadLink.onclick = function(event){
			document.body.removeChild(event.target);					
		};
		downloadLink.style.display = "none";
		document.body.appendChild(downloadLink);

		downloadLink.click();
	}
	
	compressAndSave(){
		//compress a large text string via a web worker using the lz-string library.

		console.log("Starting web worker for compression ...");
		//Make a new web worker.
		var worker = new Worker('CompressModel.js');
		var self = this;				
		worker.onmessage = function(e){
			if(e.data.compressedText){
				console.log(e.data.message+"\n---\n"+e.data.compressedText);
				self.saveToFile(e.data.compressedText);
			}
		};
		
		worker.onerror = function(e) {
			console.log('Error: Line ' + e.lineno + ' in ' + e.filename + ': ' + e.message);
		};
		
		//start the worker
		function replacer(key, value) {
		  if (typeof value === "number"){return value.round(8);}
		  return value;
		}
		var textToCompress = JSON.stringify(this.model[this.currentSolutionKey], replacer);				
		worker.postMessage({'plaintext': textToCompress,});
	}
	
	loadFromFile(){
		self = this;
		//create a new file chooser.
		//http://www.w3schools.com/jsref/dom_obj_fileupload.asp			
		var fileSelector = document.createElement('input');
		fileSelector.setAttribute('type', 'file');
		fileSelector.onchange = function(){
			//a file has been chosen, so load it.
			var reader = new FileReader();
			reader.onload = function(){
				self.decompressModel(this.result);
			};
			reader.readAsText(this.files[0]);
		};

		//simulate a click on a 'browse' button.
		fileSelector.click();
	}

	decompressModel(compressedText){
		var self = this;
		//decompress a large text string via a web worker using the lz-string library.

		//display a temporary message.
		self.updateSolverResults("<br><br><br>Loading model.<br>Please wait ...", false);

		console.log("Starting web worker for decompression ...");
		//Make a new web worker.
		var worker = new Worker('CompressModel.js');
		worker.onmessage = function(e){
			if(e.data.plaintext){
				try{
					var data = JSON.parse(e.data.plaintext);
					var x = BsplineBlade.prototype.clone(data.model);
					
					//add and set the model.
					self.add(x, data.resultsText);
					self.set(self.currentSolutionKey);
				}
				catch(err){
					console.log(err);
				}
			}
		};
		
		worker.onerror = function(e) {
			console.log('Error: Line ' + e.lineno + ' in ' + e.filename + ': ' + e.message);
		};
		
		//start the worker
		worker.postMessage({'compressedText': compressedText,});
	}

	updateSolverResults(text, hasButtons){
		this.parent.outputToScreen.solverResults.innerHTML = text;
		this.parent.outputToScreen.showSolverResults(true, hasButtons);
		if(hasButtons && this.currentSolutionKey)this.parent.refresh({gui:true, modelRedraw:true, results:true});
	}

	updateDisplay(key){
		//update display and menus.
		if(key)this.updateSolverResults(this.model[key].resultsText, true);
		else{
			this.parent.outputToScreen.showSolverResults(false);
			this.parent.refresh({gui:true, modelRedraw:true, results:false});
		}	
	}
	//end of Solutions methods
}//end of Solutions	
//
//
////
class Settings{
	constructor(x){
		this.settings = {};
		if(x)this.set(x);
	}	
	
	//methods
	add(name,x){
		if(x.value && x.defaultValue === undefined)x.defaultValue = x.value;
		this.settings[name] = x;
	}
	
	set(x){
		for(name in x)this.add(name, x[name]);
	}

	getValues(isDefault){
		var x = {};
		
		if(isDefault)this.setDefaults();
		
		for(name in this.settings){
			if(this.settings[name] && this.settings[name].value !== undefined)x[name] = this.settings[name].value;
		}	
		
		return x;
	}

	setDefaults(){
		for(name in this.settings)
			if(this.settings[name].defaultValue !== undefined)
				this.settings[name].value = this.settings[name].defaultValue;
	}
	//end of Settings methods
}//end of Settings