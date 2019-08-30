function [changedArea,pctErrL,pctErrG,nAdj,aErrTbInit] = scale_ball_volumes_divFree(areaToChange,area,ncell)
cellsInv = 1:ncell;
nCellsInv = length(cellsInv);
pctErrL = zeros(ncell,1);
aErrL = zeros(ncell,1); % L represents the error in local mass balance
scalFactorL = zeros(ncell,1);
pctErrG = zeros(ncell,1);
aErrG = zeros(ncell,1); % G represents the error in global mass balance
scalFactorG = zeros(ncell,1);
eTol = 0.05; % tolerance in pct error
maxIter = 10;
nAdj=1;
scFactorG = zeros(ncell,maxIter);
changedArea = areaToChange;
for j=1:nCellsInv
    currCellL = cellsInv(j);
    pctErrL(currCellL) = abs(sum(changedArea(:,currCellL))-area(currCellL))/area(currCellL);
    aErrL(currCellL) = sum(changedArea(:,currCellL))-area(currCellL); % if +, excess, if -, missing
end
aErrTbInit = aErrL;

for j=1:nCellsInv
    currCellG = cellsInv(j);
    pctErrG(currCellG) = abs(sum(changedArea(currCellG,:))-area(currCellG))/area(currCellG);
    aErrG(currCellG) = sum(changedArea(currCellG,:))-area(currCellG); % if +, excess, if -, missing
end

while (nAdj<=maxIter && (max(pctErrL)>eTol || max(pctErrG)>eTol))
    for j=1:nCellsInv
        currCellG = cellsInv(j);
        aErrG(currCellG) = sum(changedArea(currCellG,:))-area(currCellG); % if +, excess, if -, missing
        pctErrG(currCellG) = abs(sum(changedArea(currCellG,:))-area(currCellG))/area(currCellG);
        scalFactorG(currCellG) = area(currCellG)/sum(changedArea(currCellG,:));
    end
    scFactorG(:,nAdj) = scalFactorG;
    for i=1:nCellsInv %% cells involved/needing adjustment (global mass balance)
        currCell = cellsInv(i);
        areasInvG = find(changedArea(currCell,:)>0);
        changedArea(currCell,areasInvG) = changedArea(currCell,areasInvG)*scalFactorG(currCell);
    end
    
    for j=1:nCellsInv
        currCellL = cellsInv(j);
        pctErrL(currCellL) = abs(sum(changedArea(:,currCellL))-area(currCellL))/area(currCellL);
        scalFactorL(currCellL) = area(currCellL)/sum(changedArea(:,currCellL));
    end
    
    for i=1:nCellsInv %% cells involved/needing adjustment (local mass balance)
        currCell = cellsInv(i);
        areasInvL = find(changedArea(:,currCell)>0);
        changedArea(areasInvL,currCell) = changedArea(areasInvL,currCell)*scalFactorL(currCell);
    end
    nAdj = nAdj+1;
end
for j=1:nCellsInv
    currCellL = cellsInv(j);
    pctErrL(currCellL) = abs(sum(changedArea(:,currCellL))-area(currCellL))/area(currCellL);
    aErrL(currCellL) = sum(changedArea(:,currCellL))-area(currCellL);
end
if nAdj > maxIter
    nAdj = maxIter;
end
scFactorG = scFactorG(:,1:nAdj);
end