# spherical_gradient_paths
Gradient path generation for functions defined on a sphere

Uses the function f(t,p)=cos(p)sin(t)+cos(2t)

theta=t and phi=p is used throughout

Functions:
-----------------------------------------------
plotfunc3d() plots the toy problem function on a sphere in 3D.

plotfunc2d() plots the toy problem function in 3D in the theta-phi plane.

plotfunc2ddom() plots the toy problem function in 3D in the theta-phi plane and includes a wider domain than values that define the sphere.

plotcritpts2d() plots the toy problem function in 3D in the theta-phi plane and marks critical points.

grad(t,p) returns the gradient at the point (t,p). Uses the gradient for spherical coordinates.

gradpath(t,p) returns the gradient path that goes through the point (t,p) using Euler's method.

gradpathrk4(t,p) returns the gradient path that goes through the point (t,p) using RK4.

spheretocart(gp) converts a gradient path, gp, from spherical coordinates to Cartesian coordinates.

getphis(t) get the correct phi angles for the gradient path through the poles of the sphere by checking which value of phi gives the steepest and ascent and decent. Input values are t=0 or t=pi for the two poles. 

plotpath(t,p) plots the gradient path through the point (t,p) in 2D in the theta-phi plane using Euler's method.

plotpathrk4(t,p) plots the gradient path through the point (t,p) in 2D in the theta-phi plane using RK4.

plotpath3d(t,p) plots the gradient path through the point (t,p) in 3D on the sphere using Euler's method.

plotpath3drk4(t,p) plots the gradient path through the point (t,p) in 3D on the sphere using RK4.

gradpathflag(t,p) returns the gradient path that goes through the point (t,p) using Euler's method. Flags degenerate paths through the poles by returning a 0 if this happens.

gradpathflagrk4(t,p) returns the gradient path that goes through the point (t,p) using RK4. Flags degenerate paths through the poles by returning a 0 if this happens.

getpoints(n) generates a set of n points in the theta-phi plane to seed gradient paths.

plotallpaths2d(n) generates n paths in the theta-phi plane using Euler's method and plots all the paths that are not degenerate. Displays the amount of paths that were excluded.

plotallpaths2drk4(n) generates n paths in the theta-phi plane using RK4 and plots all the paths that are not degenerate. Displays the amount of paths that were excluded.

plotallpaths3d(n) generates n paths in 3D on the sphere using Euler's method and plots all the paths that are not degenerate. Displays the amount of paths that were excluded.

plotallpaths3drk4(n) generates n paths in 3D on the sphere using RK4 and plots all the paths that are not degenerate. Displays the amount of paths that were excluded.
