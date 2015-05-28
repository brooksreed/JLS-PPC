% Script to run a simple simulation of JLS-PPC
% BR, 5/27/2015

% SISO double integrator system
% a few options for schedules

% calls simJLSPPC (which calls functions in core)


%{
5/28
- tested with tm=0,1,2
- tested with tc=0,1,2 (checked XhMPC vs. Xh)
 
%}


clear variables
close all
clc

% (with very long ACK history, a posteriori estimate should have no effects
% of control packet losses.  control still affected due to XhMPC lacking
% information in real-time loop)
nACKHistory = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% schedule
tm = 1; % meas delay
tc = 1; % control delay
ta = 1; % ACK delay

% 'SISO_' options for _piggyback or _noACK:
% (NOTE - if piggyback, ta should equal tm)
% 'SISOALL' - [1],[1] - std. discrete time
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi 

%sched = 'SISO4_piggyback';
%sched = 'SISO4_noACK';
%sched = 'SISO2_noACK';
%sched = 'SISO2_piggyback';
sched = 'SISOALL_piggyback';

% packet success probabilities
alphaBar = .2; % controls
betaBar = 1;  % measurements
gammaBar = .8; % ACKs

%%%%

Ns = 50; % sim length
NpMult = 5; % the MPC horizon Np = Ts*NpMult 
Nv = 1;   % # vehicles (comms channels)

% Double Integrator
A =[1,1;0,1];
Bu = [0.5;1];
Bw = Bu;
C = [1 0];      % position output
%C = eye(2);    % full state output

Q = [10,0;0,1];
Qf = 10*Q;
R = 1;
umax = 10;
umin = -10;
Xmax = [];Xmin = [];
nLevels = 15;   % quantization levels
codebook = linspace(umin,umax,nLevels);

% process/measurement noise
W = .1;
V = .1;%4;
% W = 1;
% V = .1;

% estimation init
xHat1 = zeros(2,1);
% cov... uncertain position but better-known velocity (closer to zero)
P1 = [25,0;0,9];     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xIC = 5*randn(2,1);
% position only, no initial velocity (so like step resp from rest)
xIC(2)=0;xIC(1)=5;

% (IF WANT TO DEBUG CONTROLLER - INIT ESTIMATOR PERFECTLY)
% xHat1 = xIC;P1 = 1*eye(2);

w = sqrt(W)*randn(1,Ns);
v = sqrt(V)*randn(1,Ns);

% packet loss sequences
beta = zeros(Nv,Ns);alpha = zeros(Nv,Ns);gamma = zeros(Nv,Ns);
for k = 1:Ns
    beta(:,k) = (sign(rand(Nv,1) - (1-(betaBar)))*0.5 + 0.5);
    alpha(:,k) = (sign(rand(Nv,1) - (1-(alphaBar)))*0.5 + 0.5);
    gamma(:,k) = (sign(rand(Nv,1) - (1-(gammaBar)))*0.5 + 0.5);
end

[Pi,Xi,Lambda,tap,Ts] = createSchedule(sched,Nv,Ns,tc);

if(strfind(sched,'piggyback'))
    % ACK piggybacked to measurement
    gamma = beta;   % overwrite
end
Np = NpMult*Ts;

[r] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tap,...
    alphaBar,Pi,Xi,Lambda,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha,beta,gamma,nACKHistory);

% convenient if want to save:
r.P1 = P1;
r.v = v;
r.w = w;
r.xIC = xIC;
r.xHat1 = xHat1;
r.alpha = alpha;
r.beta = beta;
r.gamma = gamma;
    
    
%% plots (for SISO systems)

CPlot = [1 0];  % output for plotting
plotXhMPC = 1;  % plots prediction used for computing control
plotLosses = 1; % plots packet losses for c, m (no a right now)

% (could be used with a saved r struct too)
NxSys = size(r.P,1);    % underlying system states (no buffer)
Ns = size(r.X,2);

figure

subplot(3,1,[1 2])
hx = plot(0:Ns-1,CPlot*r.X(1:NxSys,:));
hold on
hxh = plot(0:Ns-1,CPlot*r.Xh(1:NxSys,:),'g.:');
title(sprintf('integrator sys, alphaBar = %0.2f, W=%.1f, V=%.1f',alphaBar,W,V))

hb=[];hc=[];
if(plotLosses)
    for k = 1:Ns
        if(r.alpha(k)==0)
            hc = plot([k-1+tc k-1+tc], [-0.25 0],'r');
        end
        if(r.beta(k)==0)
            hb = plot([k-1 k-1],[-.75 -.5],'m');
        end
%         if(r.gamma(k)==0)
%             % need to save/load tap, do on per-veh. basis
%             ha = plot([k-1-tc-tap k-1-tc-tap],[-.75 -.5],'c');
%         end
    end
end

if(plotXhMPC)
    hxhm = plot(0:Ns-1,CPlot*squeeze(r.XhMPC(1:NxSys,end,1:Ns)),'ko');
    legend([hx hxh hxhm hc hb],'X','XHat','XHatMPC','c loss','m loss')
else
    legend([hx hxh hc hb],'X','XHat','c loss','m loss')
end


subplot(3,1,3)
stairs(0:Ns-1,r.u,'k')
hold on
stairs(0:Ns-1,r.w,'b:')
legend('u','w')
xlabel('time step')
