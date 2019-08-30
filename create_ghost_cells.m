function [cell_n,nGhostCells,ghost_v,ghost_n,ghost_e,all_vertex,ghost_areas,ghost_centers,origGhostVert]=create_ghost_cells(ncell,nedge,nvert,vertex,diam,cell_n,cell_v,cell_e,tstep)
% creates ghost cells along the boundary of omega (to be used when tracking
% goes out of domain)
nGhostCells = 0;

for i=1:ncell % loop over all cells to locate boundary cells and add ghost cells
    nbe = length(cell_e{i}); % number of edges
    for j=1:nbe
        if cell_n{i}(j)==0
            nGhostCells = nGhostCells+1;
        end
    end
end
nGhostCells = nGhostCells+4; % square domain, 4 corners, thus 4 additional cells, need to change if domain is not square
origGhostVert = zeros(nvert,2); 
nGhostVerts = nGhostCells+4; % 1 new vertex corresponding to each new quadrilateral, +4 because 4 corners
ghost_v = cell(1,nGhostCells);
ghost_n = cell(1,nGhostCells);
ghost_e = cell(1,nGhostCells);
all_vertex = zeros(nvert+nGhostVerts,2);
all_vertex(1:nvert,:) = vertex(1:nvert,:);
g_ht = max(diam/sqrt(2))*max(tstep,1);
% g_ht = 0.0625;
LL = intersect(find(vertex(:,1)==0),find(vertex(:,2)==0));
LR = intersect(find(vertex(:,1)==1),find(vertex(:,2)==0));
UL = intersect(find(vertex(:,1)==0),find(vertex(:,2)==1));
UR = intersect(find(vertex(:,1)==1),find(vertex(:,2)==1));
allCorners = [LL LR UL UR];
ghostCellCtr = 0; % counter for labeling ghost cells
for i=1:ncell % loop over all cells to locate boundary cells and add ghost cells
    nbe = length(cell_e{i}); % number of edges
    for j=1:nbe
        if cell_n{i}(j)==0
            ghostCellCtr = ghostCellCtr+1;
            currGhostCell = ncell+ghostCellCtr; %cell no. (global) of current ghost cell
            cell_n{i}(j) = currGhostCell;
            ghost_n{ghostCellCtr} = [i 0 0 0]; %based on construction, 3rd edge is always the bdry edge
            ghost_e{ghostCellCtr} = [cell_e{i}(j) nedge+2*ghostCellCtr-1 nedge+2*ghostCellCtr 0];
            ghost_v{ghostCellCtr} = [cell_v{i}(j+1) cell_v{i}(j) nvert+ghostCellCtr 0 cell_v{i}(j+1)];
            %check if there are common vertices (will signal an
            %intersection if yes) (note: must exclude corners bec there are
            %3 ghost cells that use the corner vertices)
            %first vertex
            
            if origGhostVert(cell_v{i}(j+1),1)==0
                origGhostVert(cell_v{i}(j+1),1) = ghostCellCtr;
            else %means that one cell already used it (and the ghost cell id is stored)
                 %also means this is 2nd vertex for the other cell
                 origGhostVert(cell_v{i}(j+1),2) = ghostCellCtr;
                 if all(allCorners~=cell_v{i}(j+1))
                 ghost_n{origGhostVert(cell_v{i}(j+1),1)}(2) = currGhostCell;
                 ghost_n{ghostCellCtr}(4) = origGhostVert(cell_v{i}(j+1))+ncell;
                 ghost_e{ghostCellCtr}(4) = ghost_e{origGhostVert(cell_v{i}(j+1),1)}(2);
                 ghost_v{ghostCellCtr}(4) = ghost_v{origGhostVert(cell_v{i}(j+1),1)}(3);
                 end
            end
            
            %2nd vertex
            
            if origGhostVert(cell_v{i}(j),1)==0
                origGhostVert(cell_v{i}(j),1) = ghostCellCtr;
            else %means that one cell already used it (and the ghost cell id is stored)
                %also means this is 1st vertex for the other cell
                origGhostVert(cell_v{i}(j),2) = ghostCellCtr;
                if all(allCorners~=cell_v{i}(j))
                 ghost_n{origGhostVert(cell_v{i}(j),1)}(4) = currGhostCell;
                 ghost_e{origGhostVert(cell_v{i}(j),1)}(4) = ghost_e{ghostCellCtr}(2);
                 ghost_v{origGhostVert(cell_v{i}(j),1)}(4) = ghost_v{ghostCellCtr}(3);
                 ghost_n{ghostCellCtr}(2) = origGhostVert(cell_v{i}(j),1)+ncell;
                end
            end
            

            %identify if horizontal or vertical boundary by comparing
            %change in x and change in y
            xDiff = abs(vertex(cell_v{i}(j),1)-vertex(cell_v{i}(j+1),1));
            yDiff = abs(vertex(cell_v{i}(j),2)-vertex(cell_v{i}(j+1),2));
            %% check which boundary we are at for the square domain
            if yDiff<xDiff %horizontal bdry
                if vertex(cell_v{i}(j),1)<vertex(cell_v{i}(j+1),1) %bottom bdry
                    all_vertex(nvert+ghostCellCtr,:) = [vertex(cell_v{i}(j),1) vertex(cell_v{i}(j),2)-g_ht];
                else %top bdry
                    all_vertex(nvert+ghostCellCtr,:) = [vertex(cell_v{i}(j),1) vertex(cell_v{i}(j),2)+g_ht];
                end
            else %vertical bdry
                if vertex(cell_v{i}(j),2)<vertex(cell_v{i}(j+1),2) %right bdry
                    all_vertex(nvert+ghostCellCtr,:) = [vertex(cell_v{i}(j),1)+g_ht vertex(cell_v{i}(j),2)];
                else %left bdry
                    all_vertex(nvert+ghostCellCtr,:) = [vertex(cell_v{i}(j),1)-g_ht vertex(cell_v{i}(j),2)];
                end
            end
        end
    end
end
% create the 4 cells at the corners 
for i=1:4
    ghostCellCtr = ghostCellCtr+1;
    currCorner = allCorners(i);
    if i==1
        all_vertex(nvert+nGhostCells-4+2*i-1,:) = [vertex(currCorner,1)-g_ht vertex(currCorner,2)];
    elseif i==2
        all_vertex(nvert+nGhostCells-4+2*i-1,:) = [vertex(currCorner,1) vertex(currCorner,2)-g_ht];
    elseif i==3
        all_vertex(nvert+nGhostCells-4+2*i-1,:) = [vertex(currCorner,1) vertex(currCorner,2)+g_ht];
    elseif i==4
         all_vertex(nvert+nGhostCells-4+2*i-1,:) = [vertex(currCorner,1)+g_ht vertex(currCorner,2)];
    end
    all_vertex(nvert+nGhostCells-4+2*i,:) = [vertex(currCorner,1)+(-1)^i*g_ht vertex(currCorner,2)+((-1)^ceil(i/2))*g_ht];
    cornerLoc = find(ghost_v{origGhostVert(currCorner,1)}==currCorner); % access this to know which ghost cells are this cell's current neighbors
    ghost_v{ghostCellCtr} = [currCorner nvert+(nGhostCells-4)+2*i-1 nvert+(nGhostCells-4)+2*i 0  currCorner]; 
    ghost_e{ghostCellCtr} = [nedge+2*(nGhostCells-4)+3*i-2 nedge+2*(nGhostCells-4)+3*i-1 nedge+2*(nGhostCells-4)+3*i 0]; 
    if any(cornerLoc==1)
        ghost_v{origGhostVert(currCorner,1)}(4) = ghost_v{ghostCellCtr}(2);
        ghost_v{ghostCellCtr}(4) = ghost_v{origGhostVert(currCorner,2)}(3);
        ghost_e{ghostCellCtr}(4) = ghost_e{origGhostVert(currCorner,2)}(2);
        ghost_n{ghostCellCtr} = [origGhostVert(currCorner,1)+ncell 0 0 origGhostVert(currCorner,2)+ncell];
        ghost_n{origGhostVert(currCorner,1)}(4) = ghostCellCtr+ncell;
        ghost_e{origGhostVert(currCorner,1)}(4) = ghost_e{ghostCellCtr}(1);
        ghost_n{origGhostVert(currCorner,2)}(2) = ghostCellCtr+ncell;
    else
        ghost_v{origGhostVert(currCorner,2)}(4) = ghost_v{ghostCellCtr}(2);
        ghost_v{ghostCellCtr}(4) = ghost_v{origGhostVert(currCorner,1)}(3);
        ghost_e{ghostCellCtr}(4) = ghost_e{origGhostVert(currCorner,1)}(2);
        ghost_n{ghostCellCtr} = [origGhostVert(currCorner,2)+ncell 0 0 origGhostVert(currCorner,1)+ncell];
        ghost_n{origGhostVert(currCorner,2)}(4) = ghostCellCtr+ncell;
        ghost_e{origGhostVert(currCorner,2)}(4) = ghost_e{ghostCellCtr}(1);
        ghost_n{origGhostVert(currCorner,1)}(2) = ghostCellCtr+ncell;
    end
end
ghost_areas = zeros(nGhostCells,1);
for i=1:nGhostCells
    ghost_areas(i) = compute_area(all_vertex(ghost_v{i},:));
end
ghost_centers = gravity_centers(nGhostCells,ghost_v,all_vertex,ghost_areas);

end