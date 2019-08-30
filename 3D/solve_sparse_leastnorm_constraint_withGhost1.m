function changedArea = solve_sparse_leastnorm_constraint_withGhost1(areaToChange,area,ncell,nGhostCells);
nz = find(areaToChange>0);
changedArea = areaToChange;
numnz = (1:length(nz))';
posI = zeros(ncell*ncell,1);
posJ = zeros(ncell*ncell,1);
posZ = zeros(ncell*ncell,1);
bConst = zeros(2*ncell,1);
bConst(1:ncell) = area; % local mass conservation (need to change if non-div free velocity field)
bConst(ncell+1:2*ncell) = area; % global mass conservation
ctr=1;
for i=1:ncell
    nzCol = find(areaToChange(:,i)>0);
    nzColLoc = nzCol + (i-1)*(ncell+nGhostCells);
    for j=1:length(nzColLoc)
        %constraint for local mass cons
        if ~isempty(numnz(nz==nzColLoc(j)))
            posI(ctr)=i;
            posJ(ctr)=numnz(nz==nzColLoc(j));
            posZ(ctr)=areaToChange(nz(numnz(nz==nzColLoc(j))));
            ctr=ctr+1;
        end
    end
    nzRow = find(areaToChange(i,:)>0);
    nzRowLoc = (nzRow-1)*(ncell+nGhostCells)+i;
    for j=1:length(nzRowLoc)
        % constraint for global mass cons
        if ~isempty(numnz(nz==nzRowLoc(j)))
            posI(ctr)=i+ncell;
            posJ(ctr)=numnz(nz==nzRowLoc(j));
            posZ(ctr)=areaToChange(nz(numnz(nz==nzRowLoc(j))));
            ctr=ctr+1;
        end
    end
end
posI(posI==0)=[];
posJ(posJ==0)=[];
posZ(posZ==0)=[];
maxNo = nz(max(posJ));
nz(nz>maxNo)=[];
Aconst = sparse(posI,posJ,posZ);
H=speye(length(nz));
lb = -1+1e-4*ones(size(nz));
ub = 1-1e-4*ones(size(nz));
Aconst = Aconst(1:2*ncell-1,:);
bConst = bConst(1:2*ncell-1);
bConst = bConst -Aconst*ones(size(nz));
options = optimoptions('quadprog','Algorithm','interior-point-convex');
options.StepTolerance = 0;
xadj = quadprog(H,[],[],[],Aconst,bConst,lb,ub,[],options);

if ~isempty(xadj)
    changedArea(nz) = changedArea(nz) + xadj.*changedArea(nz);
end