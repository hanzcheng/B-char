
function [cB,r,w,mainCell] = compute_multi_radius_weights(ncell,cell_v,vertex,area,snBalls)
nBalls=snBalls^2;
bTotal=ncell*nBalls;
cB=zeros(bTotal,2); %ball center
r = zeros(bTotal,1); %radius
w = zeros(bTotal,1); %weight 
mainCell=zeros(bTotal,1);
ctr=1;
for i=1:ncell
    xL=min(vertex(cell_v{i},1));
    xR=max(vertex(cell_v{i},1));
    yL=min(vertex(cell_v{i},2));
    yR=max(vertex(cell_v{i},2));
    xDist = xR-xL; %width of the cell
    yDist = yR-yL; %length of the cell
    
    
    
    xNow=linspace(xL+(xR-xL)/(2*snBalls),xR-(xR-xL)/(2*snBalls),snBalls);
    yNow=linspace(yL+(yR-yL)/(2*snBalls),yR-(yR-yL)/(2*snBalls),snBalls);
    [x,y]=meshgrid(xNow,yNow);
    x=x(:);
    y=y(:);
    cB(ctr:ctr+nBalls-1,1)=x;
    cB(ctr:ctr+nBalls-1,2)=y;
    mainCell(ctr:ctr+nBalls-1)=i;
    r(ctr:ctr+nBalls-1) = min(xDist,yDist)/(2*snBalls); % main circle
    Acirc = pi*r(ctr)*r(ctr)*nBalls; %area of all circles inside the cell
    
    if Acirc>area(i)
        Acirc = area(i);
        r(ctr:ctr+nBalls-1) = sqrt(Acirc/(nBalls*pi));
    end
    w(ctr:ctr+nBalls-1) = area(i)/Acirc;
    ctr=ctr+nBalls;
end


