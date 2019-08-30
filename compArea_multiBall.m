function [AreaIn,currDist,tInt]=compArea_multiBall(ncell,mainCell,origCentre, newCentre, origR, newR, origW, newW)
tic;
nBalls = length(origR);
minRad = min([newR; origR]);
tol = 1e-4*minRad; %distances that are less than 0.01% of the smallest radius are considered to come from rounding errors
AreaIn=zeros(ncell,ncell); %AreaIn(i,j) gives the area of the intersection of cell i and tracked cell j
aExpect=zeros(nBalls,1);
for i=1:nBalls
    aExpect(i)=pi*newR(i)*newR(i);
    currCentre = repmat(newCentre(i,:),nBalls,1);
    currDist = (currCentre - origCentre).^2;   %compute distance bet. centres
    currDist = sqrt(sum(currDist,2));
    currDist(currDist<tol)=0;
    %           min(currDist)
    cellsInv = find(currDist<origR+newR(i)-tol); %find the balls which intersect the the tracked ball i
    %           cellsInv
    if ~isempty(cellsInv)
        if all(mainCell(cellsInv)==mainCell(cellsInv(1))) % tracked ball is fully contained in 1 residing cell, no need to compute intersection
            AreaIn(mainCell(cellsInv(1)),mainCell(i)) = AreaIn(mainCell(cellsInv(1)),mainCell(i)) + pi*newW(cellsInv(1)) * (newR(i)^2);
        else
            cellsFinN = find(currDist(cellsInv)+origR(cellsInv)<=newR(i)); % balls which are fully contained in the new ball
            cellsPin = find(currDist(cellsInv) > abs(newR(i)-origR(cellsInv))); % balls that partially intersect the new ball
            
            if ~isempty(cellsPin)
            cosC1 = 0.5*(currDist(cellsInv(cellsPin)).^2 -origR(cellsInv(cellsPin)).^2 + newR(i)^2)./(newR(i)*currDist(cellsInv(cellsPin)));
            
            anArc1 = acos(cosC1);
            aArc1 = anArc1*newR(i)^2; % area of arc (from tracked ball i)
            cosC2 = 0.5*(currDist(cellsInv(cellsPin)).^2 +origR(cellsInv(cellsPin)).^2 - newR(i)^2)./(origR(cellsInv(cellsPin)).*currDist(cellsInv(cellsPin)));
            
            anArc2 = acos(cosC2);
            
            
            aArc2 = anArc2 .* origR(cellsInv(cellsPin)).^2; % area of arc (from original ball(s) intersected by tracked ball i)
            aInt = origW(cellsInv(cellsPin)) .* (aArc1 + aArc2 - currDist(cellsInv(cellsPin))*newR(i) .* sin(anArc1)); % area(s) intersecting with tracked ball i
            aInt = aInt/sum(aInt)*aExpect(i); % distribute the expected volume of tracked ball i into these intersections proportionally
            AreaIn(mainCell(cellsInv(cellsPin)),mainCell(i)) = AreaIn(mainCell(cellsInv(cellsPin)),mainCell(i)) + aInt.*newW(cellsInv(cellsPin));
            end
            if ~isempty(cellsFinN)
                AreaIn(mainCell(cellsInv(cellsFinN)),mainCell(i)) = AreaIn(mainCell(cellsInv(cellsFinN)),mainCell(i)) + pi*newW(cellsInv(cellsFinN)) * ( newR(i)^2);
            end
        end
    end
end
tInt = toc;
end