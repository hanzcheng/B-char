% Track according to minus velocity

function xT=track_exact3D(x0,T,vel);

dt=1e-3;
Ndt=T/dt;

xT=x0;
for idt=1:Ndt
    xT = xT - dt*vel(xT(:,1),xT(:,2),xT(:,3));
%     for j=1:length(xT)
% 	xT(j,:) = xT(j,:) - dt * vel(xT(j,1),xT(j,2),xT(j,3))';
%     end
end;



