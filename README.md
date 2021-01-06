# Material Point Method
Files referenced here are located under Projects > material_point
1) Run "make" or "make all" from project root to compile the executable
	* The C++ Eigen library and Partio library are required for this project.
2) Type "./material_point" from terminal. The executable requires 1 argument--
   the name of the mesh .obj to be simulated--and has 1-3 optional arguments:
	a. You can pass in a 2nd integer argument to dictate the point density
	   of the point sampling aka the number of points to be sampled per 
	   grid cell. By default, this is 8. 
		i.e. ./material_point cube.obj 10
	b. You can pass in a 3rd integer argument to set the grid dimensions. 
	   If you pass in only these 3 arguments, the grid's x y and z dimension
	   will all be set to this 3rd argument. 
		i.e. ./material_point cube.obj 10 100
	c. You can also pass in a set of different dimensions for the grid's 
	   x, y, and z. 
		i.e. ./material_point cube.obj 10 100 64 90
   2.1) Additional parameters are located at the top of main.cpp
	- youngs_modulus
	- poissons_ratio
	- dt: simulation time step
	- offset: offset given for mesh location, useful for repositioning within grid
	- dimensions: grid dimensions (number of cells)
	- origin: minimum corner of grid
	- grid_max: maximum corner of grid
	- grid_buffer: number of cells lining the outside of the grid
	- mass: object mass
   2.2) Number of frames is set at the bottom of main.cpp; replace the number in 
	"driver.run()" with however many frames you want to simulate.
3) The executable will print out the number of particles in the simulation, based on 
   the grid dimensions and particle density set. If satisfied with the number of 
   particles sampled, typing "y" will prompt the simulation to start, and typing 
   anything else will cause the simulation to abort. 
4) Output simulation files will be in an output folder. These .bgeo files can
   be rendered in Houdini using the rendering.hipnc attached. 

The current values are optimized for rendering the following:
	./material_point data/cube.obj 15 100

This should give you a nice, bouncy jello cube! :)