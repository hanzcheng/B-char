
function [cB,r,w,mainCell] = compute_radius_weights3D(ncell,cell_v,vertex,volume,cnBalls)
nBalls=cnBalls^3;
bTotal=ncell*nBalls;
cB=zeros(bTotal,3); %ball center
r = zeros(bTotal,1); %radius
w = zeros(bTotal,1); %weight 
mainCell=zeros(bTotal,1);
ctr=1;
for i=1:ncell
    xL=min(vertex(cell_v{i},1));
    xR=max(vertex(cell_v{i},1));
    yL=min(vertex(cell_v{i},2));
    yR=max(vertex(cell_v{i},2));
    zL=min(vertex(cell_v{i},3));
    zR=max(vertex(cell_v{i},3));
    xDist = xR-xL; %width of the cell
    yDist = yR-yL; %length of the cell
    zDist = zR-zL; %height of the cell 
    
    
    xNow=linspace(xL+(xR-xL)/(2*cnBalls),xR-(xR-xL)/(2*cnBalls),cnBalls);
    yNow=linspace(yL+(yR-yL)/(2*cnBalls),yR-(yR-yL)/(2*cnBalls),cnBalls);
    zNow=linspace(zL+(zR-zL)/(2*cnBalls),zR-(zR-zL)/(2*cnBalls),cnBalls);
    [x,y,z]=meshgrid(xNow,yNow,zNow);
    x=x(:);
    y=y(:);
    z=z(:);
    cB(ctr:ctr+nBalls-1,1)=x;
    cB(ctr:ctr+nBalls-1,2)=y;
    cB(ctr:ctr+nBalls-1,3)=z;
    mainCell(ctr:ctr+nBalls-1)=i;
    r(ctr:ctr+nBalls-1) = min([xDist,yDist,zDist])/(2*cnBalls); % main circle
    Vcirc = 4/3*pi*r(ctr)*r(ctr)*r(ctr)*nBalls; %volume of all balls inside the cell
    
    if Vcirc>volume(i)
        Vcirc = volume(i);
        r(ctr:ctr+nBalls-1) = nthroot(0.75*Vcirc/(nBalls*pi),3);
    end
    w(ctr:ctr+nBalls-1) = volume(i)/Vcirc;
    ctr=ctr+nBalls;
end


