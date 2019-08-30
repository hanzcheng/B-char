function solex = compute_exactSol(v_exact,tFin,ncell,cell_v,vertex,area,center)

%% note: by choosing snQuadPts = 2^n, it is as if we have the sol. from refining the mesh by n times
snQuadPts = 4; %sqrt of the no. of quadrature points 
% [diam,area] = diameters_areas(ncell,cell_v,vertex);
[cB,~,~,mainCell] = compute_multi_radius_weights(ncell,cell_v,vertex,area,snQuadPts);
pointsToTrack = cB;
%edge midpoints

tstep=1e-1;
nbSteps=ceil(tFin/tstep);

solini=initCond(center);
solex = zeros(size(solini));

% write_solution_ensight(0,0,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);
% write_solution_ensight(1,solex,0,times,ncell,nedge,nvert,cell_v,cell_n,cell_e,vertex,cg);

solex_ref = initCond(track_exact(pointsToTrack,nbSteps*tstep,v_exact)); % "exact" sol. on a refined grid
solex_ref(solex_ref==1)=1;
for i=1:ncell
    cellsNow = find(mainCell==i);
    nCellsNow = length(cellsNow);
    solex(i) = sum(solex_ref(cellsNow))/nCellsNow;
end

end