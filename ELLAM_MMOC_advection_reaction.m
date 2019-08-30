%% Characteristic-based schemes for the advection-reaction problem
% for solenoidal fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_t + alpha*c + div (cV) = f
% alpha is the reaction coefficient
% V is the velocity with div(V)=0
% with source term f and initial condition c(x,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
%
% -------------------------------------------------------------------
% Setup and load mesh
mesh = {'mesh2_3.mat'};       % The mesh file in MATLAB data format
loadmesh=strcat('load ../../matlab_meshes/',mesh{1});
eval(loadmesh);
cg = gravity_centers(ncell,cell_v,vertex,area);
% -------------------------------------------------------------------
% Setup problem parameters
% reaction coefficient
alfa = 0;
% Source term
source = @(x,y) 0;
global testcase
testcase = 4;

% Exact velocity and its inflow/outflow at boundary
if testcase <= 2
    v_exact = @(x,y) [ 0.0625 *ones(size(x))  0 * ones(size(y))];
    noFlow = 0; % some inflow or outflow at boundary
elseif testcase >= 3 
    v_exact = @(x,y) [(1-2*y) .* x .* (1-x), -(1-2*x) .* y .* (1-y)];
    noFlow = 1;  %no inflow or outflow at boundary
end

% Initial condition

cPrev = zeros(ncell+nedge,1);
for i=1:ncell
    cPrev(i) = initCond(cg(i,:));
end

%% choice of numerical scheme and time step
numScheme = "ELLAM"; %actually since we look at solenoidal fields
% numScheme = "MMOC"; %ELLAM and MMOC are equivalent
tstep = 0.8; %time step
tFin = 8; %final time
solex = compute_exactSol(v_exact,tFin,ncell,cell_v,vertex,area,center); % benchmark / "exact" solution, for 0 reaction and source term
nbSteps = ceil(tFin/tstep);
times=[0:tstep:tFin];
snBalls = 2;
if noFlow==0 %if there is inflow/outflow, create ghost cells at the boundary before proceeding
    [cell_n,nGhostCells,ghost_v,ghost_n,ghost_e,all_vertex,ghost_areas,ghost_centers] = create_ghost_cells(ncell,nedge,nvert,vertex,diam,cell_n,cell_v,cell_e,tstep);
    allCells = 1:(ncell+nGhostCells);
    all_v = [cell_v ghost_v];
    all_e = [cell_e ghost_e];
    all_n = [cell_n ghost_n];
    all_areas = [area;ghost_areas];
    [cB,origR,origW,mainCell] = compute_multi_radius_weights(ncell+nGhostCells,all_v,all_vertex,all_areas,snBalls);
else %otherwise, compute the radius and weights/densities of each ball
    [cB,origR,origW,mainCell] = compute_multi_radius_weights(ncell,cell_v,vertex,area,snBalls);
end
newR = origR; %solenoidal field, so volumes should not change after transport
newW = origW; %solenoidal field, so the amt of mass inside the cell should not change



pointsToTrack = cB;

% characteristic tracking
if numScheme == "ELLAM"
    [newPoints,tTrack]=complete_tracking_exact_vel(tstep,pointsToTrack,v_exact,0);
elseif numScheme == "MMOC"
    [newPoints,tTrack]=complete_tracking_exact_vel(tstep,pointsToTrack,v_exact,1);
end
% area of intersection bet tracked and residing cells

if noFlow==0 % with ghost cells
    [AreaIn,currDist,tInt]=compArea_multiBall(ncell+nGhostCells,mainCell,pointsToTrack, newPoints, origR, newR, origW, newW);
else % without ghost cells
    [AreaIn,currDist,tInt]=compArea_multiBall(ncell,mainCell,pointsToTrack, newPoints, origR, newR, origW, newW);
end

tic;
% volume adj for local and global mass conservation
AreaIn(abs(AreaIn)<1e-10)=0;
[AreaIn,pctErr,errG,nAdj,aErrTbInit] = scale_ball_volumes_divFree(AreaIn,area,ncell);

if nAdj>1
    if noFlow==0 % with ghost cells
        AreaIn = solve_sparse_leastnorm_constraint_withGhost1(AreaIn,area,ncell,nGhostCells);
    else % without ghost cells
        AreaIn = solve_sparse_leastnorm_constraint1(AreaIn,area,ncell);
    end
end
tAdj = toc;
tTotal = tAdj+tInt+tTrack;
tTrack = tTrack/tTotal;
tAdj = tAdj/tTotal;
tInt = tInt/tTotal;
write_solution_ensight(0,0,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
%process the initial condition for plotting by
% setting the value at the edges to be the average of the value of 2
% adjacent cells
for i=1:ncell
    nbe = length(cell_e{i});
    for j=1:nbe
        if cPrev(ncell+cell_e{i}(j))==0
            cPrev(ncell+cell_e{i}(j)) = cPrev(i);
        else
            cPrev(ncell+cell_e{i}(j)) =0.5*( cPrev(ncell+cell_e{i}(j))+cPrev(i));
        end
    end
end

write_solution_ensight(1,cPrev,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);

%measure of mass conservation
massCons = zeros(1,nbSteps);
%% solve the advection PDE
for m=1:nbSteps
    A = zeros(ncell,1);
    b = zeros(ncell,1);
    
    % assemble LHS and source term
    for cell_i=1:ncell
        A(cell_i) = area(cell_i);
        if numScheme == "ELLAM"
            intC = AreaIn(1:ncell,cell_i)'*cPrev(1:ncell);
        elseif numScheme == "MMOC"
            intC = AreaIn(cell_i,1:ncell)*cPrev(1:ncell);
        end
        b(cell_i) = intC +tstep* (source(cg(cell_i,1),cg(cell_i,2)) - alfa*cPrev(cell_i));
    end
    
    % solve for c
    cCurr = b ./ A;
    
    %measure of global mass conservation bet. current time and prev. time
    massCons(m) = sum(area.*cPrev(1:ncell))-sum(area.*cCurr(1:ncell));
    
    cPrev = cCurr;
    c = zeros(ncell+nedge,1);
    c(1:ncell) = cPrev;
    % set the value at the edges to be the average of the value of 2
    % adjacent cells (to be used for plotting)
    for i=1:ncell
        nbe = length(cell_e{i});
        for j=1:nbe
            if c(ncell+cell_e{i}(j))==0
                c(ncell+cell_e{i}(j)) = cCurr(i);
            else
                c(ncell+cell_e{i}(j)) =0.5*( c(ncell+cell_e{i}(j))+cCurr(i));
            end
        end
    end
    c(abs(c)<1e-12)=0;
    
    write_solution_ensight(1,c,m,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
end
[L2err,L1err,pL2err,pL1err,allLayers,nLayers] = compute_errors(solex,cCurr,area,ncell);
L2err
L1err
pL2err
pL1err

if numScheme == "ELLAM"
    write_solution_vtk(c,'solution_ELLAM_rect',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
elseif numScheme == "MMOC"
    write_solution_vtk(c,'solution_MMOC_rect',ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex);
end
