function trackedPoints = FO_Euler_exactVel3D(pointsToTrack,dtTrack,dv_exact,tracking)
%% Performs characteristic tracking by using first order euler
% inputs: hho_FE -- the structure of the finite element
%       : pointsToTrack -- the set of points that we track
%       : currCell(i) -- points to the cell whose velocity will be used to 
%                       track pointsToTrack(i,:)
%       : dtTrack -- time step
%       : RT_k -- the RT_k FE used to approximate the velocity field
%       : tracking -- informs us whether we forward or backward track (0 is
%       backward tracking)
% output: trackedPoints -- the location of pointsToTrack after tracking

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tracking ==0
    dtTrack=-dtTrack;
end
x = pointsToTrack(:,1);
y = pointsToTrack(:,2);
z = pointsToTrack(:,3);

k_1 = dv_exact(x,y,z);
k_1 = k_1 .*[dtTrack dtTrack dtTrack];
x1 = x+k_1(:,1);
y1 = y+k_1(:,2);
z1 = z+k_1(:,3);
trackedPoints = [x1 y1 z1];
% trackedPoints(abs(trackedPoints(:,1)-1000)<bdTol,1)=1000;
% trackedPoints(abs(trackedPoints(:,1))<bdTol,1)=0;
% trackedPoints(abs(trackedPoints(:,2)-1000)<bdTol,2)=1000;
% trackedPoints(abs(trackedPoints(:,2))<bdTol,2)=0;
