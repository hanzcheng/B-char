ball-approximated ELLAM or MMOC scheme for an advection-reaction problem:

Main file: 

ELLAM_MMOC_advection_reaction - solves an advection-reaction problem (with given velocity) via an ELLAM or MMOC scheme

Initialisation: 
1. gravity_centers - computes the cell center of gravity
2. compute_multi_radius_weights- gives the center(s) of each ball, together with its radius
3. initCond - gives the initial conditions used for the test
5. compute_exactSol - computes a benchmark "exact" solution by using method of characteristics and very small time steps
5. track_exact - tracking used for the benchmark "exact" solution
6. create_ghost_cells - creates a set of ghost cells at the boundary (to be used when we have inflow / outflow on the boundary)


Files related to characteristic tracking: 

1. FO_Euler_exactVel - performs the characteristic tracking by using a forward Euler scheme and using micro timesteps
2. complete_tracking_exact_vel - performs the full tracking (i.e. loops and perform FO_Euler_exactVel until the time step is exhausted)
3. scale_ball_volumes_divFree - scales the volumes of the balls to reduce the error in global mass conservation
4. solve_sparse_leastnorm_constraint - solves the optimisation problem: mimimise the volume change subject to the constraint of local and global mass conservation
(takes advantage of using sparse representation, so that the computations are quick)
(_withghost is used when there is inflow/outflow at the boundary) 

Processing:
1. compute_area - computes the area of a given cell
2. compArea_multiBall - computes area of intersection bet. balls
3. compute_errors - computes the errors in L^1 and L^2, and also gives the number of layers for which the numerical solution and benchmark solution differ by more than 5%
4. write_solution_vtk - writes a vtk file for visualising the solution profile
5. write_solution_ensight - writes an ensight file for visualising the solution profile (note: Here, the results are written in a different folder, in this case, specified to be ensight/Un. Users may opt to change where the results are written by changing lines 121-134.) 

A sample test result is given in the folder ensight. This is obtained by running ELLAM_MMOC_advection_reaction as is. Of course, the mesh, time step, source term, velocity field, and initial conditions may be modified to run different test cases on different types of meshes. Convergence can also be tested by refining the mesh and reducing the time step.

Note: the functions have the same description in 3D. However, in 3D, the mesh is limited to Cartesian type meshes.