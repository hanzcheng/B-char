function solex = compute_exactSol3D(v_exact,tFin,ncell,cell_v,vertex,volume,center,testCase)

%% note: by choosing snQuadPts = 2^n, it is as if we have the sol. from refining the mesh by n times
cnQuadPts = 4; %cbrt of the no. of quadrature points 
% [diam,area] = diameters_areas(ncell,cell_v,vertex);
[cB,~,~,mainCell] = compute_radius_weights3D(ncell,cell_v,vertex,volume,cnQuadPts);
pointsToTrack = cB;
%edge midpoints

tstep=1e-1;
nbSteps=ceil(tFin/tstep);

solini=initCond3D(center,testCase);
solex = zeros(size(solini));

% write_solution_ensight(0,0,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
% write_solution_ensight(1,solex,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);

solex_ref = initCond3D(track_exact3D(pointsToTrack,nbSteps*tstep,v_exact),testCase); % "exact" sol. on a refined grid
solex_ref(solex_ref==1)=1;
for i=1:ncell
    cellsNow = find(mainCell==i);
    nCellsNow = length(cellsNow);
    solex(i) = sum(solex_ref(cellsNow))/nCellsNow;
end

end