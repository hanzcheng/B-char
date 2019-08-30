% create a 3D grid with cubes
function [ncell,nface,nvert,face_v,cell_f,cell_n,cell_v,vertex,volume,center]=create_3D_grid(nCubes)

ncell = nCubes^3;%no. of cells/cubes
sLen = 1/nCubes; %side length of each cube
x=0:sLen:1;
y=0:sLen:1;
z=0:sLen:1;
[X,Y,Z]=meshgrid(x,y,z);
xxx=X(:); % x, y, and z coordinates of the vertices
yyy=Y(:);
zzz=Z(:);
vertex=[xxx yyy zzz]; %vertices
nvert = length(vertex); %no. of vertices
cell_f = cell(1,ncell);
cell_v = cell(1,ncell); %vertices on a given cell
cell_n = cell(1,ncell); %neighbors of a given cell
volume = zeros(ncell,1);
center = zeros(ncell,3); %cell barycenters
volume = volume + sLen^3; %since made up of all cubes of equal lengths, this is the volume (otherwise can compute one X x Y x Z through a loop)
nface = (3*(nCubes+1))*nCubes^2;
%% alternative way of computing no. of faces
% nface1 = 0 ;
% while cubeCtr>=2
%     if cubeCtr==nCubes
%         nface1 = nface1 + 20*cubeCtr-20;
%     else
%         nface1 = nface1 + 16*cubeCtr-20;
%     end
%     cubeCtr = cubeCtr-2;
% end
% if cubeCtr == 1
%     if cubeCtr == nCubes
%         nface1 = 6;
%     else
%         nface1 = nface1+2;
%     end
% end
% nface1 = nCubes*nface1 - nCubes*nCubes*(nCubes-1);
face_v = cell(nface,1); %vertices on a given face'
fPerLayer = nCubes^2; % no. of faces in each layer (also same as no. of cells in each layer)
nLayers = nCubes+1;
currBdry = 0;
fTol = sLen*1e-4; %tolerance used when searching for points
for j=1:nLayers
    faceCtr = 1;
    for i=1:fPerLayer %loop over all the faces in each layer
        currLayerX = find(abs(vertex(:,1)-currBdry)<fTol);
        currLayerY = find(abs(vertex(:,2)-currBdry)<fTol);
        currLayerZ = find(abs(vertex(:,3)-currBdry)<fTol);
        currLevel = floor(faceCtr/nCubes);
        if mod(faceCtr,nCubes)==0
            currLevel = currLevel-1;
        end
        %cutting through z (fill the base before moving to next level)
        
        maxXz = sLen*(i-currLevel*nCubes);
        if maxXz == 0
            maxXz = 1;
        end
        maxYz = sLen*(currLevel+1);
        if maxYz>1
            maxYz = 1;
        end
        Z_vert1 = intersect(find(abs(vertex(:,1)-maxXz)<fTol),find(abs(vertex(:,2)-maxYz)<fTol));
        Z_vert1 = intersect(Z_vert1,currLayerZ);
        Z_vert2 = intersect(find(abs(vertex(:,1)-maxXz)<fTol),find(abs(vertex(:,2)-(maxYz-sLen))<fTol));
        Z_vert2 = intersect(Z_vert2,currLayerZ);
        Z_vert3 = intersect(find(abs(vertex(:,1)-(maxXz-sLen))<fTol),find(abs(vertex(:,2)-(maxYz-sLen))<fTol));
        Z_vert3 = intersect(Z_vert3,currLayerZ);
        Z_vert4 = intersect(find(abs(vertex(:,1)-(maxXz-sLen))<fTol),find(abs(vertex(:,2)-maxYz)<fTol));
        Z_vert4 = intersect(Z_vert4,currLayerZ); % note: counterclockwise orientation when viewed as pointing towards the negative direction
        face_v{i+(j-1)*fPerLayer} = [Z_vert1 Z_vert2 Z_vert3 Z_vert4 Z_vert1];
        %cutting through y (fill the front before going to the next level)
        maxXy = sLen*(i-currLevel*nCubes);
        if maxXy == 0
            maxXy = 1;
        end
        maxZy = sLen*(currLevel+1);
        if maxZy >1
            maxZy = 1;
        end
        Y_vert1 = intersect(find(abs(vertex(:,1)-maxXy)<fTol),find(abs(vertex(:,3)-maxZy)<fTol));
        Y_vert1 = intersect(Y_vert1,currLayerY);
        Y_vert2 = intersect(find(abs(vertex(:,1)-(maxXy-sLen))<fTol),find(abs(vertex(:,3)-maxZy)<fTol));
        Y_vert2 = intersect(Y_vert2,currLayerY);
        Y_vert3 = intersect(find(abs(vertex(:,1)-(maxXy-sLen))<fTol),find(abs(vertex(:,3)-(maxZy-sLen))<fTol));
        Y_vert3 = intersect(Y_vert3,currLayerY);
        Y_vert4 = intersect(find(abs(vertex(:,1)- maxXy)<fTol),find(abs(vertex(:,3)-(maxZy-sLen))<fTol));
        Y_vert4 = intersect(Y_vert4,currLayerY); % note: counterclockwise orientation when viewed as pointing towards the negative direction
        face_v{i+(j-1)*fPerLayer+nLayers*fPerLayer} = [Y_vert1 Y_vert2 Y_vert3 Y_vert4 Y_vert1];
        %cutting through x (fill the left before moving to next level)
        maxYx = sLen*(i-currLevel*nCubes);
        if maxYx == 0
            maxYx = 1;
        end
        maxZx = sLen*(currLevel+1);
        if maxZx >1
            maxZx = 1;
        end
        X_vert1 = intersect(find(abs(vertex(:,2)-maxYx)<fTol),find(abs(vertex(:,3)-maxZx)<fTol));
        X_vert1 = intersect(X_vert1,currLayerX);
        X_vert2 = intersect(find(abs(vertex(:,2)-(maxYx-sLen))<fTol),find(abs(vertex(:,3)-maxZx)<fTol));
        X_vert2 = intersect(X_vert2,currLayerX);
        X_vert3 = intersect(find(abs(vertex(:,2)-(maxYx-sLen))<fTol),find(abs(vertex(:,3)-(maxZx-sLen))<fTol));
        X_vert3 = intersect(X_vert3,currLayerX);
        X_vert4 = intersect(find(abs(vertex(:,2)-maxYx)<fTol),find(abs(vertex(:,3)-(maxZx-sLen))<fTol));
        X_vert4 = intersect(X_vert4,currLayerX); % note: counterclockwise orientation when viewed as pointing towards the negative direction
        face_v{i+(j-1)*fPerLayer+2*nLayers*fPerLayer} = [X_vert1 X_vert2 X_vert3 X_vert4 X_vert1];
        faceCtr = faceCtr + 1;
    end
    currBdry = currBdry + sLen;
end
%cell faces, in the ff order: bottom, top, front, back, left, right
currCell = 1;
cellsInv = zeros(nface,1); % cells that interact with a given face
for i=1:nCubes %move on to next level (i.e. from bottom to top)
    for j=1:nCubes %front to back
        for k=1:nCubes % fill in base from left to right
            cell_n{currCell} = zeros(1,6);
            cell_bot = k+(j-1)*nCubes+(i-1)*fPerLayer;
            cell_top = k+(j-1)*nCubes+i*fPerLayer;
            if cellsInv(cell_bot)==0
                cellsInv(cell_bot) = currCell;
            else
                cell_n{currCell}(1) = cellsInv(cell_bot);
                cell_n{cellsInv(cell_bot)}(2) = currCell;
            end
            cellsInv(cell_top) = currCell;
            cell_front = k+(j-1)*fPerLayer + (i-1)*nCubes + nLayers*fPerLayer;
            cell_back = k+(j-1)*fPerLayer + (i-1)*nCubes + (nLayers+1)*fPerLayer;
            if cellsInv(cell_front)==0
                cellsInv(cell_front) = currCell;
            else
                cell_n{currCell}(3) = cellsInv(cell_front);
                cell_n{cellsInv(cell_front)}(4) = currCell;
            end
            cellsInv(cell_back) = currCell;
            cell_left = j+(k-1)*fPerLayer+(i-1)*nCubes+2*nLayers*fPerLayer;
            cell_right = j+(k-1)*fPerLayer+(i-1)*nCubes+(2*nLayers+1)*fPerLayer;
            if cellsInv(cell_left)==0
                cellsInv(cell_left) = currCell;
            else
                cell_n{currCell}(5) = cellsInv(cell_left);
                cell_n{cellsInv(cell_left)}(6) = currCell;
            end
            cellsInv(cell_right) = currCell;
            cell_f{currCell} = [cell_bot cell_top cell_front cell_back cell_left cell_right];
            currCell = currCell+1;
        end
    end
end
for i=1:ncell
    cell_v{i} = [face_v{cell_f{i}(1)} face_v{cell_f{i}(2)}]; %vertices are given from bottom face, then top face
    center(i,:) = 0.125*(sum(vertex(face_v{cell_f{i}(1)}(1:4),:))+sum(vertex(face_v{cell_f{i}(2)}(1:4),:)));
end