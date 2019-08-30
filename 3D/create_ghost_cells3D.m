function [cell_n,nGhostCells,ghost_v,ghost_n,ghost_f,all_vertex,ghost_volumes]=create_ghost_cells3D(ncell,nvert,vertex,sLen,cell_n,face_v,cell_f,nCubes,tstep)

% creates ghost cells along the boundary of omega (to be used when tracking
% goes out of domain)
nGhostCells = 6*nCubes^2; % only create ghost cells on direct neighbors
fTol = sLen*1e-4;
sAdj = sLen*max(tstep,1);
nGhostVerts = 6*(nCubes+1)^2; % creating an nxn grid on each boundary face 
ghost_v = cell(1,nGhostCells);
ghost_n = cell(1,nGhostCells);
ghost_f = cell(1,nGhostCells);
all_vertex = zeros(nvert+nGhostVerts,3);
all_vertex(1:nvert,:) = vertex(1:nvert,:);
all_vertex(nvert+1:nvert+(nCubes+1)^2,:) = vertex(abs(vertex(:,3))<fTol,:)-repmat([0 0 sAdj],(nCubes+1)^2,1); %bot
all_vertex(nvert+1*(nCubes+1)^2+1:nvert+2*(nCubes+1)^2,:) = vertex(abs(vertex(:,3)-1)<fTol,:) + repmat([0 0 sAdj],(nCubes+1)^2,1); %top
all_vertex(nvert+2*(nCubes+1)^2+1:nvert+3*(nCubes+1)^2,:) = vertex(abs(vertex(:,2))<fTol,:) - repmat([0 sAdj 0],(nCubes+1)^2,1); %front
all_vertex(nvert+3*(nCubes+1)^2+1:nvert+4*(nCubes+1)^2,:) = vertex(abs(vertex(:,2)-1)<fTol,:) + repmat([0 sAdj 0],(nCubes+1)^2,1); %back
all_vertex(nvert+4*(nCubes+1)^2+1:nvert+5*(nCubes+1)^2,:) = vertex(abs(vertex(:,1))<fTol,:) - repmat([sAdj 0 0],(nCubes+1)^2,1); %left
all_vertex(nvert+5*(nCubes+1)^2+1:nvert+6*(nCubes+1)^2,:) = vertex(abs(vertex(:,1)-1)<fTol,:) + repmat([sAdj 0 0],(nCubes+1)^2,1); %right

ghostCtr = 0;
for i=1:ncell
    nbf = size(cell_f{i},2);
    for j=1:nbf
    if cell_n{i}(j)==0
        ghostCtr = ghostCtr+1;
        currVerts = zeros(1,5);
        if j==1
            allZ = find(abs(all_vertex(:,3)+sAdj)<fTol);
            currVertY = vertex(face_v{cell_f{i}(j)},2); %y coord of current vertex
            currVertX = vertex(face_v{cell_f{i}(j)},1);
            for k=1:4
            vertLoc = intersect(allZ,find(abs(all_vertex(:,2)-currVertY(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,1)-currVertX(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [currVerts face_v{cell_f{i}(j)}];
        elseif j==2
            allZ = find(abs(all_vertex(:,3)-1-sAdj)<fTol);
            currVertY = vertex(face_v{cell_f{i}(j)},2); %y coord of current vertex
            currVertX = vertex(face_v{cell_f{i}(j)},1);
            for k=1:4
            vertLoc = intersect(allZ,find(abs(all_vertex(:,2)-currVertY(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,1)-currVertX(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [face_v{cell_f{i}(j)} currVerts];
        elseif j==3
            allY = find(abs(all_vertex(:,2)+sAdj)<fTol);
            currVertZ = vertex(face_v{cell_f{i}(j)},3); %y coord of current vertex
            currVertX = vertex(face_v{cell_f{i}(j)},1);
            for k=1:4
            vertLoc = intersect(allY,find(abs(all_vertex(:,3)-currVertZ(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,1)-currVertX(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [currVerts face_v{cell_f{i}(j)}];
        elseif j==4 %back neighbor
            allY = find(abs(all_vertex(:,2)-1-sAdj)<fTol);
            currVertZ = vertex(face_v{cell_f{i}(j)},3); %y coord of current vertex
            currVertX = vertex(face_v{cell_f{i}(j)},1);
            for k=1:4
            vertLoc = intersect(allY,find(abs(all_vertex(:,3)-currVertZ(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,1)-currVertX(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [face_v{cell_f{i}(j)} currVerts];
        elseif j==5 %left
            allX = find(abs(all_vertex(:,1)+sAdj)<fTol);
            currVertZ = vertex(face_v{cell_f{i}(j)},3); %y coord of current vertex
            currVertY = vertex(face_v{cell_f{i}(j)},2);
            for k=1:4
            vertLoc = intersect(allX,find(abs(all_vertex(:,3)-currVertZ(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,2)-currVertY(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [currVerts face_v{cell_f{i}(j)}];
        else
            allX = find(abs(all_vertex(:,1)-1-sAdj)<fTol);
            currVertZ = vertex(face_v{cell_f{i}(j)},3); %y coord of current vertex
            currVertY = vertex(face_v{cell_f{i}(j)},2);
            for k=1:4
            vertLoc = intersect(allX,find(abs(all_vertex(:,3)-currVertZ(k))<fTol));
            vertLoc = intersect(vertLoc,find(abs(all_vertex(:,2)-currVertY(k))<fTol));
            currVerts(k) = vertLoc;
            end
            currVerts(5) = currVerts(1);
            ghost_v{ghostCtr} = [currVerts face_v{cell_f{i}(j)}];
        end
    end
    end
end
ghost_volumes = zeros(nGhostCells,1);
for i=1:nGhostCells
    ghost_volumes(i) = sAdj^3; %since we are doing a domain of unif. cubes
end

end