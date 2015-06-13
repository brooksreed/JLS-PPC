% Script to run a simple simulation of JLS-PPC
% BR, 5/27/2015

% double integrator system
% a few options for schedules
% handles varying lengths of ACK histories

% calls simJLSPPC (which calls functions in core)
% some simple plots for SISO systems at end

clear variables
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sys = 'DoubleIntegrator';
sys = 'SCALAR';

% schedule

% 'SISO_' options for _piggyback or _noACK:
% (NOTE - if piggyback, ta should equal tm)
% 'SISOALL' - [1],[1] - std. discrete time
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi 

%sched = 'SISO4_piggyback';
%sched = 'SISO4_noACK';
%sched = 'SISO2_noACK';
%sched = 'SISO2_piggyback';

%sched = 'SISO2ALLCONTROL_noACK';
%sched = 'SISO2ALLCONTROL_piggyback';
sched = 'SISOALL_piggyback';
%sched = 'SISOALL_noACK';

% # ACK Histories sent (makes most sense to be multiple of schedule length)
% For 'SINGLE ACK': set nACKHistory = Ts (schedule length)
nACKHistory = 5;

% NOTE:
% With very long ACK history, a posteriori estimate should have no effects
% of control packet losses.  Control is still affected due to XhMPC 
% lacking information in the real-time loop

% delays
tm = 1; % meas delay
tc = 1; % control delay
ta = 1; % ACK delay

% packet success probabilities
alpha_cBar = .5; % controls
alpha_mBar = .7;  % measurements
alpha_aBar = .7; % ACKs (if piggyback used, betaBar overrides gammaBar)
covPriorAdj = 1;

%%%%%%%%%%%%%%%%%%

Ns = 20; % sim length
NpMult = 4; % the MPC horizon Np = Ts*NpMult 
Nv = 1;   % # vehicles (comms channels)

switch sys
    case 'DoubleIntegrator'
        
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
        W = 1;
        V = 1;%4;
        % W = 1;
        % V = .1;
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P1 = [25,0;0,9];     

    case 'SCALAR'
        
        % Scalar integrator
        A =1;
        Bu = 1;
        Bw = Bu;
        C = 1;      % position output

        Q = 10;
        Qf = 10*Q;
        R = 1;
        umax = 1;
        umin = -1;
        Xmax = [];Xmin = [];
        nLevels = 33;   % quantization levels
        codebook = linspace(umin,umax,nLevels);

        % process/measurement noise
        W = 1;
        V = 1;%4;
        % W = 1;
        % V = .1;
        
        P1 = 9;
        
end

% estimation init
xHat1 = zeros(size(A,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xIC = 5*randn(size(A,1),1);

if(strfind(sys,'DoubleIntegrator'))
    % position only, no initial velocity (so like step resp from rest)
    xIC(2)=0;xIC(1)=5;
end

% (IF WANT TO DEBUG CONTROLLER - INIT ESTIMATOR PERFECTLY)
% xHat1 = xIC;P1 = 1*eye(2);

w = sqrt(W)*randn(size(Bw,1),Ns);
v = sqrt(V)*randn(size(C,1),Ns);

% packet loss sequences
alpha_m = zeros(Nv,Ns);alpha_c = zeros(Nv,Ns);alpha_a = zeros(Nv,Ns);
for k = 1:Ns
    alpha_m(:,k) = (sign(rand(Nv,1) - (1-(alpha_mBar)))*0.5 + 0.5);
    alpha_c(:,k) = (sign(rand(Nv,1) - (1-(alpha_cBar)))*0.5 + 0.5);
    alpha_a(:,k) = (sign(rand(Nv,1) - (1-(alpha_aBar)))*0.5 + 0.5);
end

% hardcoded sequences for consistent debugging
%alpha_c(:,1:11) = [1 1 0 0 1 1 0 1 0 0 1];
alpha_m(:,1:11) =  [1 1 1 0 1 1 0 0 0 1 1];

[Pi_c,Pi_m,Pi_a,tac,Ts] = createSchedule(sched,Nv,Ns,tc);

if(strfind(sched,'piggyback'))
    % ACK piggybacked to measurement
    alpha_a = alpha_m;   % overwrite
end
Np = NpMult*Ts;

%% call sim fcn

[r] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tac,...
    alpha_cBar,Pi_c,Pi_m,Pi_a,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha_c,alpha_m,alpha_a,covPriorAdj,nACKHistory);

% convenient if want to save:
r.P1 = P1;
r.v = v;
r.w = w;
r.xIC = xIC;
r.xHat1 = xHat1;
r.alpha_c = alpha_c;
r.alpha_m = alpha_m;
r.alpha_a = alpha_a;
% (could reconstruct from sched but simpler to just save)
r.Pi_c = Pi_c;
r.Pi_m = Pi_m;
r.Pi_a = Pi_a;
    
% (save r struct here if want)

%% plots (for SISO systems)

if(size(A,1)==2)
    CPlot = [1 0];  % output for plotting
else
    CPlot = 1;
end

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
title(sprintf('integrator sys, alphaBar = %0.2f, W=%.1f, V=%.1f',alpha_cBar,W,V))

hb=[];hc=[];
if(plotLosses)
    for k = 1:Ns
        if( (r.Pi_c(k) - r.alpha_c(k)*r.Pi_c(k)) == 1 )
            %hc = plot([k-1+tc k-1+tc], [-0.25 0],'r');
            hc = plot(k-1+tc,0,'ro');
        end
        if( (r.Pi_m(k) - r.alpha_m(k)*r.Pi_m(k)) ==1 )
            %hb = plot([k-1 k-1],[-.75 -.5],'m');
            hb = plot(k-1,-2,'mo');
        end
%         if( (r.Pi_a(k) - r.alpha_a(k)*r.Pi_a(k)) ==1 )
%             % need to save/load tac, do on per-veh. basis
%             ha = plot([k-1-tc-tac k-1-tc-tac],[-.75 -.5],'c');
%         end
    end
end

if(plotXhMPC)
    hxhm = plot(0:Ns-1,CPlot*squeeze(r.XhMPC(1:NxSys,end,1:Ns)),'k.');
    legend([hx hxh hxhm hc hb],'X','XHat','XHatMPC','c loss','m loss')
else
    legend([hx hxh hc hb],'X','XHat','c loss','m loss')
end

subplot(3,1,3)
stairs(0:Ns-1,r.u,'k')
hold on
stairs(0:Ns-1,r.uNoLoss,'b:')
legend('u true','u planned')
%stairs(0:Ns-1,r.w,'b:')
%legend('u','w')
xlabel('time step')

if(size(r.P,1)==1)
    figure
    plot(squeeze(r.P))
end

