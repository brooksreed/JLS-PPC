% sets up JLSPPC system for simulations
% pulls these settings outside of run script

% BR, 6/18/2015

%%%%%%%%%%%%%%%%%%

if(strfind(sched,'piggyback'))
    if(ta~=tm)
        disp('warning - resetting ta to tm')
        ta = tm;
    end
    alpha_aBar = alpha_mBar;
elseif(strfind(sched,'NoACK'))
    ta = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM INIT/SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch system
    case 'SISO_DOUBLE_INTEGRATOR'
        
        % controller settings
        umax = 10;      
        umin = -10;
        nLevels = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W = 1;
        V = 1;
                
        % cov... uncertain position but better-known velocity (closer to zero)
        P1 = [25,0;0,9];     
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = [0.5;1];
        Bw = Bu;
        C = [1 0];      % position output
        Xmax = [];Xmin = [];
        codebook = linspace(umin,umax,nLevels);
        
 case 'MIMO_DOUBLE_INTEGRATOR'
     
        % controller settings
        umax = [10;10];
        umin = [-10;-10];
        nLevels = 15;   % quantization levels
        Q = [10,0;0,1];
        Qf = 10*Q;
        R = eye(2);
               
        % process/measurement noise
        W = eye(2);
        V = eye(2);
        
        % cov... uncertain position but better-known velocity (closer to zero)
        P1 = [25,0;0,9];  
        
        % Double Integrator
        A =[1,1;0,1];
        Bu = eye(2);
        Bw = Bu;
        C = eye(2);    % full state output
        Xmax = [];Xmin = [];
        codebook = linspace(min(umin),max(umax),nLevels);
        
    case 'SCALAR'
        
        % controller settings        
        umax = 1;
        umin = -1;
        nLevels = 33;   % quantization levels
        Q = 10;
        Qf = 10*Q;
        R = 1;
        
        % process/measurement noise
        W = 1;
        V = 1;
        
        % initial covariance
        P1 = 9;
        
        % Scalar integrator
        A =1;
        Bu = 1;
        Bw = Bu;
        C = 1;      % position output
        Xmax = [];Xmin = [];
        codebook = linspace(umin,umax,nLevels);
        
end

% check
if( length(alpha_cBar)~=Nv )
    disp('WARNING: alpha_c, Nv mismatch')
end
if( length(alpha_mBar)~=Nv )
    disp('WARNING: alpha_m, Nv mismatch')
end
if( rem(size(C,1),Nv) )
    disp('WARNING: Nv not into C')
end
if( rem(size(Bu,2),Nv) )
    disp('WARNING: Nv not into Bu')
end    
    
% estimation init
xHat1 = zeros(size(A,1),1);

% random time-series realizations
w = sqrt(W)*randn(size(Bw,2),Ns);
v = sqrt(V)*randn(size(C,1),Ns);

% packet loss sequences
alpha_m = zeros(Nv,Ns);alpha_c = zeros(Nv,Ns);alpha_a = zeros(Nv,Ns);
for k = 1:Ns
    alpha_m(:,k) = (sign(rand(Nv,1) - (1-(alpha_mBar)))*0.5 + 0.5);
    alpha_c(:,k) = (sign(rand(Nv,1) - (1-(alpha_cBar)))*0.5 + 0.5);
    alpha_a(:,k) = (sign(rand(Nv,1) - (1-(alpha_aBar)))*0.5 + 0.5);
end

% schedule time series
[Pi_c,Pi_m,Pi_a,tac,Ts] = createSchedule(sched,Nv,Ns,tc);
Np = NpMult*Ts;

