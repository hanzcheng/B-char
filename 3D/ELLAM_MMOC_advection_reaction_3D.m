%% Characteristic-based schemes for the advection-reaction problem
% for solenoidal fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_t + alpha*c + div (cV) = f
% alpha is the reaction coefficient
% V is the velocity with div(V)=0
% with source term f and initial condition c(x,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
%
% -------------------------------------------------------------------
% Setup problem parameters
% reaction coefficient
alfa = 0;
% Source term
source = @(x,y,z) 0;
testCase = 3;
% Exact velocity and its inflow/outflow at boundary
if testCase<=2
    v_exact = @(x,y,z) [ 0.0625*ones(size(x))  0*ones(size(y))   0*ones(size(z)) ];
elseif testCase>=3
    v_exact = @(x,y,z) [(1-2*y) .* x .* (1-x), -(1-2*x) .* y .* (1-y), 0.0625*ones(size(z))];
end
noFlow = 0;
% noFlow = 1;
% final time
tFin = 8;
% -------------------------------------------------------------------
% numerical scheme and time step
numScheme = "ELLAM";
% numScheme = "MMOC";
tstep = 0.8;
nbsteps = ceil(tFin/tstep);
% -------------------------------------------------------------------
% Setup and load mesh
nCubes = 16;
sLen = 1/nCubes;
[ncell,nface,nvert,face_v,cell_f,cell_n,cell_v,vertex,volume,center]=create_3D_grid(nCubes); %create a nCubes x nCubes x nCubes grid
% initial condition and benchmark solution

solex = compute_exactSol3D(v_exact,tFin,ncell,cell_v,vertex,volume,center,testCase);  % benchmark / "exact" solution, for 0 reaction and source term
cIni = zeros(ncell,1);
for i=1:ncell
    cIni(i) = initCond3D(center(i,:),testCase);
end
cPrev = cIni;
% load the balls inscribed in each cell
cnBalls = 2; % cnBalls is cbrt of number of balls, so we have cnBalls^3 balls per cell

if noFlow ==0 %if there is inflow/outflow, create ghost cells at the boundary before proceeding
    [cell_n,nGhostCells,ghost_v,ghost_n,ghost_f,all_vertex,ghost_volumes]=create_ghost_cells3D(ncell,nvert,vertex,sLen,cell_n,face_v,cell_f,nCubes,tstep);
    allCells = 1:(ncell+nGhostCells);
    all_v = [cell_v ghost_v];
    all_f = [cell_f ghost_f];
    all_n = [cell_n ghost_n];
    all_volumes = [volume;ghost_volumes];
    [cB,origR,origW,mainCell] = compute_radius_weights3D(ncell+nGhostCells,all_v,all_vertex,all_volumes,cnBalls);
else   %otherwise, compute the radius and weights/densities of each ball
    [cB,origR,origW,mainCell] = compute_radius_weights3D(ncell,cell_v,vertex,volume,cnBalls);
end
newR = origR; %solenoidal field, so volumes should not change after transport
newW = origW; %solenoidal field, so the amt of mass inside the cell should not change



pointsToTrack = cB;

% characteristic tracking
if numScheme == "ELLAM"
    [newPoints,tTrack]=complete_tracking_exactVel3D(tstep,pointsToTrack,v_exact,0);
elseif numScheme == "MMOC"
    [newPoints,tTrack]=complete_tracking_exactVel3D(tstep,pointsToTrack,v_exact,1);
end

if noFlow==0 % with ghost cells
    [VolIn,currDist,tInt]=compArea_multiBall3D(ncell+nGhostCells,mainCell,pointsToTrack, newPoints, origR, newR, origW, newW);
else % without ghost cells
    [VolIn,currDist,tInt]=compArea_multiBall3D(ncell,mainCell,pointsToTrack, newPoints, origR, newR, origW, newW);
end

tic;
% volume adj for local and global mass conservation
VolIn(abs(VolIn)<1e-10)=0;
[VolIn,pctErr,scalFactor,nAdj,aErrTbInit] = scale_ball_volumes_divFree(VolIn,volume,ncell);
%
if nAdj>1
    if noFlow==0 % with ghost cells
        VolIn = solve_sparse_leastnorm_constraint_withGhost1(VolIn,volume,ncell,nGhostCells);
    else % without ghost cells
        VolIn = solve_sparse_leastnorm_constraint1(VolIn,volume,ncell);
    end
end
tAdj = toc;
tTotal = tAdj+tInt+tTrack;
% tTrack = tTrack/tTotal;
% tAdj = tAdj/tTotal;
% tInt = tInt/tTotal;

massCons = zeros(1,nbsteps);
%% solve the advection PDE
for m=1:nbsteps
    A = zeros(ncell,1);
    b = zeros(ncell,1);
    
    % assemble LHS and source term
    for cell_i=1:ncell
        A(cell_i) = volume(cell_i);
        if numScheme == "ELLAM"
            intC = VolIn(1:ncell,cell_i)'*cPrev(1:ncell);
        elseif numScheme == "MMOC"
            intC = VolIn(cell_i,1:ncell)*cPrev(1:ncell);
        end
        b(cell_i) = intC +tstep* (source(center(cell_i,1),center(cell_i,2),center(cell_i,3)) - alfa*cPrev(cell_i));
    end
    
    % solve for c
    cCurr = b ./ A;
    
    %measure of global mass conservation bet. current time and prev. time
    massCons(m) = sum(volume.*cPrev(1:ncell))-sum(volume.*cCurr(1:ncell));
    
    cPrev = cCurr;
end

[L2err,L1err,pL2err,pL1err,allLayers,nLayers] = compute_errors(solex,cCurr,volume,ncell);
L2err
L1err
pL2err
pL1err