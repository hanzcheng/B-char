function [VolIn,currDist,tInt]=compArea_multiBall3D(ncell,mainCell,origCentre, newCentre, origR, newR, origW, newW)
tic;
nBalls = length(origR);
minRad = min([newR; origR]);
tol = 1e-4*minRad; %distances that are less than 0.01% of the smallest radius are considered to come from rounding errors
VolIn=zeros(ncell,ncell); %VolIn(i,j) gives the volume of the intersection of cell i and tracked cell j
vExpect=zeros(nBalls,1);
for i=1:nBalls
    vExpect(i)=4/3*pi*newR(i)*newR(i)*newR(i);
    currCentre = repmat(newCentre(i,:),nBalls,1);
    currDist = (currCentre - origCentre).^2;   %compute distance bet. centres
    currDist = sqrt(sum(currDist,2));
    currDist(currDist<tol)=0;
    %           min(currDist)
    cellsInv = find(currDist<origR+newR(i)-tol); %find the balls which intersect the the tracked ball i
    %           cellsInv
    if ~isempty(cellsInv)
        if all(mainCell(cellsInv)==mainCell(cellsInv(1))) % tracked ball is fully contained in 1 residing cell, no need to compute intersection
            VolIn(mainCell(cellsInv(1)),mainCell(i)) = VolIn(mainCell(cellsInv(1)),mainCell(i)) + 4/3* pi*newW(cellsInv(1)) * (newR(i)^3);
        else
            cellsFinN = find(currDist(cellsInv)+origR(cellsInv)<=newR(i)); % balls which are fully contained in the new ball
            cellsPin = find(currDist(cellsInv) > abs(newR(i)-origR(cellsInv))); % balls that partially intersect the new ball
            
            if ~isempty(cellsPin) %volume of cap1 + volume of cap 2
            h2 = (newR(i)+origR(cellsInv(cellsPin))-currDist(cellsInv(cellsPin))).*(newR(i)-origR(cellsInv(cellsPin))+currDist(cellsInv(cellsPin))) ./ (2*currDist(cellsInv(cellsPin)));
            h1 = (newR(i)+origR(cellsInv(cellsPin))-currDist(cellsInv(cellsPin))).*(origR(cellsInv(cellsPin))-newR(i)+currDist(cellsInv(cellsPin))) ./ (2*currDist(cellsInv(cellsPin)));
            vCap1 = pi/3* h1.^2 .*(3*newR(i)-h1);
            vCap2 = pi/3* h2.^2 .*(3*origR(cellsInv(cellsPin))-h2);
            vInt = (vCap1+vCap2) .* origW(cellsInv(cellsPin)); % volume(s) intersecting with tracked ball i
            vInt = vInt/sum(vInt)*vExpect(i); % distribute the expected volume of tracked ball i into these intersections proportionally
            VolIn(mainCell(cellsInv(cellsPin)),mainCell(i)) = VolIn(mainCell(cellsInv(cellsPin)),mainCell(i)) + vInt.*newW(cellsInv(cellsPin));
            end
            if ~isempty(cellsFinN)
                VolIn(mainCell(cellsInv(cellsFinN)),mainCell(i)) = VolIn(mainCell(cellsInv(cellsFinN)),mainCell(i)) + 4/3*pi*newW(cellsInv(cellsFinN)) * (newR(i)^3);
            end
        end
    end
end
tInt = toc;
end