function [qt,zPhi,errq] = pursuitExpMeas(qx,qy,pb,gb,pt)
% converts vehicle position in local xy to scalar q and zPhi
% function [qt,zPhi,errq] = pursuitExpMeas(qx,qy,pb,gb,pt)
% qx, qy: scalar x and y vehicle position in 2D local coords
% pb, gb: (2x1) reference position and gradient
% pt: scalar ptilde perturbation (true)
%
% projects (qx,qy) onto C2 line to give qt (along C2) and errq
% makes true zPhi  based on qt and pt (true instance perturbation)
% currently always uses gradients of 1 (could easily be modified)

% BR, 9/16/2014

% angle thetaG from +x axis, CCW (RH coord sys) to the gradient line
% -pi:pi
thetaG = atan2(gb(2),gb(1));
if(thetaG<0)
    thetaG = 2*pi + thetaG;
end

% q errors relative to pb origin
dqx = qx - pb(1);
dqy = qy - pb(2);

% now rotate (qx,qy) about pb by -thetaG to put in the "perturb frame"
ang = -thetaG;
Rot = [cos(ang),-sin(ang);sin(ang),cos(ang)];
tmp = Rot*[dqx;dqy];
qt = tmp(1);
errq = tmp(2);

% make zPhi
G = 1;
zPhi = G*(pt-qt);
