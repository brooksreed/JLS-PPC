% Script to run a simple simulation of JLS-PPC
% BR, 5/27/2015

% SISO double integrator system
% a few options for schedules

% calls simJLSPPC (which calls functions in core)


clear variables
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% schedule
tm = 1; % meas delay
tc = 0; % control delay
ta = 0; % ACK delay

% 'SISO_' options for _piggyback or _noACK:
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi 

%sched = 'SISO4_piggyback';
%sched = 'SISO4_noACK';
sched = 'SISO2_noACK';

% packet success probabilities
alphaBar = 1; % controls
betaBar = 1;  % measurements
gammaBar = 1; % ACKs

%%%%

Ns=100; % sim length
NpMult = 4; % the MPC horizon Np = Ts*NpMult 
Nv=1;   % # vehicles (comms channels)

% SISO Double Integrator
A =[1,1;0,1];
Bu = [0.5;1];
Bw = Bu;
C = [1 0];
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
V = 4;
% W = 1;
% V = .1;

% estimation init
xHat1 = zeros(2,1);
P1 = 30*eye(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xIC = randn(2,1);
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

[results] = simJLSPPC(Ns,Np,A,Bu,Bw,C,Q,Qf,R,W,V,tm,tc,ta,tap,...
    alphaBar,Pi,Xi,Lambda,umax,umin,codebook,Xmax,Xmin,xIC,P1,xHat1,...
    w,v,alpha,beta,gamma);

