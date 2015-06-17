% Script to run a simple simulation of JLS-PPC
% BR, 5/27/2015

% v1.0 6/13/2015

% TO DO: 
% clean up/separate inputs vs. other? 
% one toggle for 'debugging mode'?
% better plots
% automated saving/logging?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double integrator or scalar systems
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

system = 'SISO_DOUBLE_INTEGRATOR';
%system = 'SCALAR';

% schedule

% 'SISO_' options for _piggyback or _noACK:
% (NOTE - if piggyback, ta should equal tm)
% 'SISOALL' - [1],[1] - std. discrete time
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi 

%sched = 'SISO4_piggyback';
%sched = 'SISO4_noACK';

%sched = 'SISO2_noACK';
sched = 'SISO2_piggyback';

%sched = 'SISO2ALLCONTROL_noACK';
%sched = 'SISO2ALLCONTROL_piggyback';
%sched = 'SISOALL_piggyback';
%sched = 'SISOALL_noACK';

% # ACK Histories sent (makes most sense to be multiple of schedule length)
% For 'SINGLE ACK': sys nACKHistory = Ts (schedule length)
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
alpha_cBar = .75; % controls
alpha_mBar = .7;  % measurements
alpha_aBar = .7; % ACKs (if piggyback used, betaBar overrides gammaBar)
covPriorAdj = 1;

%%%%%%%%%%%%%%%%%%

Ns = 80; % sim length
NpMult = 4; % the MPC horizon Np = Ts*NpMult 
Nv = 1;   % # vehicles (comms channels)

switch system
    case 'SISO_DOUBLE_INTEGRATOR'
        
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

if(size(A,1)==2)
    % position only, no initial velocity (so like step resp from rest)
    xIC(2)=0;xIC(1)=5;
end

% (IF WANT TO DEBUG CONTROLLER - INIT ESTIMATOR PERFECTLY)
% xHat1 = xIC;P1 = 1*eye(2);

w = sqrt(W)*randn(size(Bw,2),Ns);
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
%alpha_m(:,1:11) =  [1 1 1 0 1 1 0 0 0 1 1];

[Pi_c,Pi_m,Pi_a,tac,Ts] = createSchedule(sched,Nv,Ns,tc);

if(strfind(sched,'piggyback'))
    % ACK piggybacked to measurement
    alpha_a = alpha_m;   % overwrite
end
Np = NpMult*Ts;

%% call sim fcn

% hack to overwrite covPriorAdj and rerun with same everything else
covPriorAdj = 0;

[r] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tac,...
    alpha_cBar,Pi_c,Pi_m,Pi_a,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha_c,alpha_m,alpha_a,covPriorAdj,nACKHistory);

% add to r struct for saving:
% run-specific parameters
r.P1 = P1;
r.v = v;
r.w = w;
r.xIC = xIC;
r.xHat1 = xHat1;
r.alpha_c = alpha_c;
r.alpha_m = alpha_m;
r.alpha_a = alpha_a;
r.Pi_c = Pi_c;
r.Pi_m = Pi_m;
r.Pi_a = Pi_a;

% system
sys.system = system;
sys.sched = sched;
sys.A = A;
sys.Bu = Bu;
sys.C = C;
sys.umax = umax;
sys.umin = umin;
sys.nLevels = nLevels;
sys.Q = Q;
sys.Qf = Qf;
sys.R = R;
sys.NpMult = NpMult;
sys.W = W;
sys.V = V;
sys.alpha_cBar = alpha_cBar;
sys.alpha_mBar = alpha_mBar;
sys.alpha_aBar = alpha_aBar;
sys.ta = ta;
sys.tc = tc;
sys.tm = tm;

r.sys = sys;

% (save r struct here if want)

%% plots
plotJLSPPC_SISO(r)




