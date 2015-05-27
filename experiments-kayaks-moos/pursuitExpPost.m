function [x2D,xh2D,q2D] = pursuitExpPost(pb,gb,pt,ptHat,qt)
% converts pt and phHat perturbs into 2D x,y (for posting to viewer)
% function [x2D,xh2D,q2D] = pursuitExpPost(pb,gb,pt,ptHat);
% pb: pBar reference for this point and time (2 x 1)
% gb: gBar reference for this point and time (2 x 1)
% pt: scalar ptilde for this point and time
% ptHat: scalar estimated ptildeHat for this point and time
% qt: scalar (projected) vehicle position along C2 line

% BR, 9/16/2014

thetaG = atan2(gb(2),gb(1));
if(thetaG<0)
    thetaG = 2*pi + thetaG;
end

% now rotate (qx,qy) about pb by -thetaG to put in the "perturb frame"
ang = thetaG;
Rot = [cos(ang),-sin(ang);sin(ang),cos(ang)];

x2D = pb + Rot*[pt;0];
xh2D = pb + Rot*[ptHat;0];
q2D = pb + Rot*[qt;0];
