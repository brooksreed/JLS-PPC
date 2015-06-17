% Script to run a simple simulation of JLS-PPC
% BR, 5/27/2015

% v1.0 6/13/2015
% v1.1 6/16/2015

% TO DO: 

% clean up/separate inputs vs. other? 
% one toggle for 'debugging mode'?
% more detailed plots for SISO

% MIMO: cleaner system setup, make Nc and Nm vs. Nv? 
% MIMO basic plots

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

%system = 'SISO_DOUBLE_INTEGRATOR';
system = 'MIMO_DOUBLE_INTEGRATOR';
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
%sched = 'SISO2_piggyback';

%sched = 'SISO2ALLCONTROL_noACK';
%sched = 'SISO2ALLCONTROL_piggyback';
%sched = 'SISOALL_piggyback';
%sched = 'SISOALL_noACK';

% 'MIMO' options: MX, IL
sched = 'MX_piggyback';
%sched = 'MX_noACK';

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
if( ~isempty(strfind(system,'SISO')) || ~isempty(strfind(system,'SCALAR')))
    alpha_cBar = .75; % controls
    alpha_mBar = .7;  % measurements
    alpha_aBar = .7; % ACKs (if piggyback used, betaBar overrides gammaBar)
else
    alpha_cBar = [.75;.75];
    alpha_mBar = [.75;.75];
    alpha_aBar = [.75;.75];
end
covPriorAdj = 1;

%%%%%%%%%%%%%%%%%%

Ns = 20; % sim length
NpMult = 4; % the MPC horizon Np = Ts*NpMult 
Nv = size(alpha_cBar,1);   % # vehicles (comms channels)

switch system
    case 'SISO_DOUBLE_INTEGRATOR'
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = [0.5;1];
        Bw = Bu;
        C = [1 0];      % position output
        umax = 10;      
        umin = -10;
        Xmax = [];Xmin = [];
        nLevels = 15;   % quantization levels
        codebook = linspace(umin,umax,nLevels);
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = 1;

        % process/measurement noise
        W = 1;
        V = 1;%4;
        % W = 1;
        % V = .1;
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P1 = [25,0;0,9];     
        
 case 'MIMO_DOUBLE_INTEGRATOR'
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = eye(2);
        Bw = Bu;
        C = eye(2);    % full state output

        Q = [10,0;0,1];
        Qf = 10*Q;
        R = eye(2);
        umax = [10;10];
        umin = [-10;-10];
        Xmax = [];Xmin = [];
        nLevels = 15;   % quantization levels
        codebook = linspace(min(umin),max(umax),nLevels);

        
        % process/measurement noise
        W = eye(2);
        V = eye(2);
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P1 = [25,0;0,9];     
        covPriorAdj = 0;    % always (for now...)
        
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

% (pull-out run-specific parameters for re-running?)

[r] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tac,...
    alpha_cBar,Pi_c,Pi_m,Pi_a,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha_c,alpha_m,alpha_a,covPriorAdj,nACKHistory);

r.sys.sched = sched;
r.sys.system =system;
%r.sys.alpha_cBar = alpha_cBar; (added inside)
r.sys.alpha_mBar = alpha_mBar;
r.sys.alpha_aBar = alpha_aBar;


% (save r struct here if want)
%{

old = cd('C:\Brooks\Dropbox\Research Dropbox\MATLAB Code\JLS-PPC local');
fname = sprintf('results_%s',dateString('DHM'));
save(fname,'r')

%or rename r 
uniquer = r;
save(fname,'uniquer')
cd(old)

%}


%% plots
if(size(C,1)==1 && size(Bu,2)==1)
    plotJLSPPC_SISO(r)
end




