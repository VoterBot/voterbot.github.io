# BladeSolver
## B-spline Wind Turbine Blade with a First-Order Panel Method Aerodynamic Solver.

Bi-cubic B-spline representation scheme for wind turbine blades. The default blade is an approximation to chord and twist distributions first published in: ***Performance and wake measurements on a 3 m diameter horizontal axis wind turbine: comparison of theory, wind tunnel and field test data.* MB Anderson, DJ Milborrow, NJ Ross - 1982.**
	
Black lines are a truncated polynomial approximation to the Anderson et. al. blade design, which used a standard Blade Element Method to calculate an optimal chord/twist distribution.  The default design shown here used a low-speed SD7062 aerofoil and was developed by Professor David H. Wood at the University of Newcastle's department of mechanical engineering, and was manufactured by his company Aerogenesis Australia Pty. Ltd.
	
Red lines are a B-spline least-squares curve fit to those distributions, with the green line being the curve's control polygon.  You can drag the control points (green rectangles) to change the default distributions.
	
Blue lines are the resulting cross-sectional distributions from the NURBS blade model.  The closer the fit to the black line, the closer the NURBS model will be to the chord and twist distributions of Professor Wood and Anderson et. al.

The solver is a first-order Panel Method, a boundary element collocation method that divides the continuous surface of the blade into discrete quadrilateral surface panels, then solves the Laplace equation for velocity potential at centre of each panel.  The influence of the trailing helical wake is included, with the dipole singularity on each wake panel being equal to the circulation around the blade strip which shed it.  Central differences of the velocity potential gives the velocity field at blade collocation points, thence surface pressure, surface forces, torque and power.

You can use the menu to discretise the blade and wake in various ways.  Fewer panels makes for significantly faster solution times.  The Chrome Javascript engine seems to be much faster than the one used by Firefox.

NOTE:  To run locally, Chrome requires this switch on the command line: 

--allow-file-access-from-files

Otherwise, it gets angry at you trying to load data from a local iframe.  Firefox will bork too, but there's no switch.
Should work fine when served by a web server in any browser.
