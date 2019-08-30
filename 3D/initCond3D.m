function uv=initCond3D(pt,testcase)
x=pt(:,1); 
y=pt(:,2);
z=pt(:,3);
uv=zeros(size(pt,1),1);
if  testcase == 1 % cubic block (0.0625,0.0625,0.0625)x(0.3125,0.3125,0.3125)
      nZx1=find(x>=0.0625);
    nZx2=find(x<=0.3125);
    nZx=intersect(nZx1,nZx2);
    nZy1=find(y>=0.0625);
    nZy2=find(y<=0.3125);
    nZy=intersect(nZy1,nZy2);
    nZ=intersect(nZx,nZy);
    nZz1=find(z>=0.0625);
    nZz2=find(z<=0.3125);
    nZz=intersect(nZz1,nZz2);
    nZ=intersect(nZ,nZz);
    
uv(nZ)=1;
    
elseif testcase == 2

    uv = exp(-10*((x-.1875) .^ 2 + (y-.1875) .^ 2 + (z-.1875) .^ 2));
    uv(x<0) = 0;
    uv(y<0) = 0;
    uv(z<0) = 0;
elseif testcase == 3 %cylinder
 nZxy = find((x-.25).^2+(y-.75).^2 < 1/64);
    nZz1=find(z>=0.0625);
    nZz2=find(z<=0.3125);
    nZz=intersect(nZz1,nZz2);
    nZ=intersect(nZxy,nZz);
    uv(nZ)=1;
  
elseif testcase == 4 
    
    uv = exp(-10*((x-.25) .^ 2 + (y-.75) .^ 2+ (z-.1875) .^ 2));
    uv(x<0) = 0;
    uv(y<0) = 0;
    uv(z<0) = 0;
    uv(z>1) = 0;
end

