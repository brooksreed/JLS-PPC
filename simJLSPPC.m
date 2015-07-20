function [results] = simJLSPPC(SIM_LEN,N_HORIZON,A_SYS,Bu_SYS,Bw_SYS,...
    C_SYS,QMPC,QfMPC,RMPC,W_KF,V_KF,TAU_M,TAU_C,TAU_A,TAU_AC,ALPHAC_BAR,...
    PI_C,PI_M,PI_A,T_S,U_MAX,U_MIN,CODEBOOK,X_MAX,X_MIN,X_IC,P_1,X_HAT1,...
    w_t,v_t,alpha_c,alpha_m,alpha_a,cov_prior_adj,N_ACKHISTORY,print_debug)
% runs simulation of JLSPPC
% [results] = simJLSPPC(SIM_LEN,N_HORIZON,A_SYS,Bu_SYS,Bw_SYS,C_SYS,...
%    QMPC,QfMPC,RMPC,W_KF,V_KF,TAU_M,TAU_C,TAU_A,TAU_AC,ALPHAC_BAR,...
%    PI_C,PI_M,PI_A,T_S,U_MAX,U_MIN,CODEBOOK,XMAX,XMIN,X_IC,P_1,X_HAT1,...
%    w_t,v_t,alpha_c,alpha_m,alpha_a,cov_prior_adj,N_ACKHISTORY,print_debug)
%
% inputs approximately match variables in paper
% cov_prior_adj = {1,0} toggle for using cov. prior adjustments in KF
% print_debug = {1,0} 'global' toggle for debug printouts
% results: large struct of results

% BR, 4/23/2014

% TO DO:
% Make sure prep for KF is ready for MIMO and multiple P*s
% Add logging of P*, P**, etc.?
% Add more tests/checks?

% more verbose debug printouts from KF
print_debug_KF = 0;

% log Pstar
logPstar = 1;

% INITIALIZATION
NX_SYS = size(A_SYS,1);
N_VEH = length(ALPHAC_BAR);
N_CONTROLS_ALL = size(Bu_SYS,2);    % all controls
N_CONTROLS_VEH = N_CONTROLS_ALL/N_VEH;         % per-vehicle
N_W = size(Bw_SYS,2);
N_Y_ALL = size(C_SYS,1);
N_Y_VEH = N_Y_ALL/N_VEH;

% initialize tNoACK 
t_NoACK = zeros(N_VEH,SIM_LEN);
% first period of controls are known
[~,tmp] = find(PI_C(:,1:T_S));
tmp = tmp+TAU_AC+(TAU_C+TAU_A);
for i = 1:N_VEH
    t_NoACK(i,tmp(i):end) = 1:(SIM_LEN-tmp(i)+1);
end

% after first planned control RX, unknown
% initialize based on 2nd pi_c + tau_c ?? (need tau_ac too?)
% tmp = find(pi_c);
% start = tmp(2)+tc
% tNoACK(:,start:Ns) = 1:(Ns-start);
    %repmat(1:SIM_LEN,[N_VEH,1]);

% tNoACK *AS KNOWN EACH STEP*
t_NoACK_save = cell(1,SIM_LEN);

Pstar_save = cell(1,SIM_LEN);
if(N_VEH==1);Pstar_overwrite = zeros(1,SIM_LEN);end

if(cov_prior_adj)
    disp('COV PRIOR ADJUST ON')
    % compute Pstar terms for specific alpha_cBar
    if(length(ALPHAC_BAR)==1)
        
        % load in the symbolic lookup table - faster
        [Pstar_coefficients,Pstar_final_coefficients] = ...
            evaluatePstars(ALPHAC_BAR);
        cov.Pstar_coefficients = Pstar_coefficients;
        cov.Pstar_final_coefficients = Pstar_final_coefficients;
        
        % optionally, if nStars>11 desired, can compute numerically for a
        % specific alpha_cBar (may take a while)
        %{
        nStars = 12;
        PstarCoefficients = computePstars(nStars,alpha_cBar,'NUMERIC');
        %}
        
    else
        %
        disp('USING MULTIVEHICLE ONE-STEP P* APPROX')
    end
    
end

% turn control constraints into time series
% (could be run with time-varying if desired)
U_MAX_T = repmat(U_MAX,1,SIM_LEN+N_HORIZON);
U_MIN_T = repmat(U_MIN,1,SIM_LEN+N_HORIZON);

% Nearest-Neighbor quantizer: make partition
if(~ischar(CODEBOOK))
    PARTITION = ( CODEBOOK(2:end)+CODEBOOK(1:(end-1)) )/2;
end

% various sytem inits/preallocations
M = makeM(N_CONTROLS_VEH,N_HORIZON,N_VEH);
etmp = [1,zeros(1,N_HORIZON-1)];
E = kron(eye(N_VEH),kron(etmp,eye(N_CONTROLS_VEH)));
BW = [Bw_SYS;zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,N_W)];

X = zeros(NX_SYS+N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN);     % includes x and b
y = zeros(N_Y_ALL,SIM_LEN);                             % true y 
yh = zeros(N_Y_ALL,SIM_LEN);                            % y into estimator 
Xh = NaN*zeros(NX_SYS+N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN);   % includes xHat and bHat
P = zeros(NX_SYS,NX_SYS,SIM_LEN);
u = zeros(N_CONTROLS_ALL,SIM_LEN);
U = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN);

% initialize jump variable matrices
Dc = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN);    
Dm = zeros(N_VEH*N_Y_VEH,N_VEH*N_Y_VEH,SIM_LEN);
Da = zeros(N_VEH,N_VEH,SIM_LEN);

% estimated alpha_c, Dc, no loss controls/buffers
Dc_hat = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN+TAU_C);   
Dc_no_loss = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN+TAU_C);   
alphac_hat = repmat(ALPHAC_BAR,[1 SIM_LEN]);
for t = (TAU_C+1):SIM_LEN
    Dc_hat(:,:,t) = makeDc(PI_C(:,t-TAU_C),alphac_hat(:,t-TAU_C),N_CONTROLS_VEH,N_HORIZON);
    Dc_no_loss(:,:,t) = makeDc(PI_C(:,t-TAU_C),PI_C(:,t-TAU_C),N_CONTROLS_VEH,N_HORIZON);
end
u_no_loss = zeros(N_CONTROLS_ALL,SIM_LEN);
b_no_loss = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,SIM_LEN);

% includes xHatMPC and bHatMPC
XhMPC = NaN*zeros(N_VEH*N_HORIZON*N_CONTROLS_VEH+NX_SYS,TAU_M+TAU_C+1,SIM_LEN);	

% variables for saving MPC output/timing
Jcomp = NaN*zeros(1,SIM_LEN);
X_plan = NaN*zeros(NX_SYS,N_HORIZON,SIM_LEN);
MPC_time = NaN*zeros(SIM_LEN,1);
loop_time = zeros(SIM_LEN,1);
MPC_fail = NaN*zeros(SIM_LEN,1);

% initial state x_1 and z_1
X(1:NX_SYS,1) = X_IC;
y(:,1) = C_SYS*X_IC+v_t(:,1);

% initial buffer - all zeros
X(NX_SYS+1:end,1) = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,1);

% first step propagation - gives x_2
X(1:NX_SYS,2) =  A_SYS*X_IC + Bw_SYS*w_t(:,1);
y(:,2) = C_SYS*X(1:NX_SYS,1) + v_t(:,2);
u(:,1) = E*X(NX_SYS+1:end,1);

% Loop starts at "estimator" at time t=tm+1, when first meas (t=1) is RX'd
for t = (TAU_M+1):(SIM_LEN-1)
    looptic = tic;
    
    % determine which measurements are available at estimator at this step
    % at step t, meas. are sent at t-tm
    Dm(:,:,t-TAU_M) = makeDm(PI_M(:,t-TAU_M),alpha_m(:,t-TAU_M),N_Y_VEH);
    yh(:,t-TAU_M) = Dm(:,:,t-TAU_M)*y(:,t-TAU_M);     
    
    if(print_debug)
        fprintf('\n~~~STEP t=%d AT ESTIMATOR~~~\n',t)
        for i = 1:N_VEH
            if(PI_M(i,t-TAU_M)*alpha_m(i,t-TAU_M)==1)
                fprintf('\nt=%d, Meas %d (sent at %d) RX success\n',t,i,t-TAU_M)
            end
        end
    end  
    
    % determine ACKs available at this step
    % update Dh (and alphaHat), KFstart
    % uses tNoACK(t-ta-1) to determine how far back to update
    % updates tNoACK(t-ta) and history based on ACKs RX'd now
    % also increments a lookahead of tNoACK(t-ta+1 --> future)
    [Dc_hat,alphac_hat,Da,KF_start,t_NoACK,~] = JLSJumpEstimator(Dc_hat,...
        PI_C,Da,alpha_c,alphac_hat,PI_A,T_S,alpha_a,t,TAU_M,TAU_C,...
        TAU_A,TAU_AC,N_CONTROLS_VEH,N_HORIZON,t_NoACK,N_ACKHISTORY,print_debug);
    t_NoACK_save{t} = t_NoACK;
    
    if(print_debug)
        disp_start = t-8;
        disp_end = t+T_S+2;
        if(disp_start<1);disp_start=1;end
        if(disp_end>SIM_LEN);disp_end=SIM_LEN;end
        fprintf('\n tNoACK(%d:%d):\n',disp_start,disp_end)
        disp(t_NoACK(:,disp_start:disp_end))
    end
    
    %%%%%%%%%%%%%%%
    % run estimator
    %%%%%%%%%%%%%%%
    
    % compute buffer and control action if all control packets through
    if( (t-TAU_M-1) >= 1)
        if(t-TAU_M-1>1)
            b_prev = b_no_loss(:,t-TAU_M-2);
        else
            b_prev = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,1);
        end
        b_no_loss(:,t-TAU_M-1) = M*(eye(N_HORIZON*N_CONTROLS_VEH*N_VEH)-Dc_no_loss...
            (:,:,t-TAU_M-1))*b_prev + Dc_no_loss(:,:,t-TAU_M-1)*...
            U(:,t-TAU_M-1);
        u_no_loss(:,t-TAU_M-1) = E*b_no_loss(:,t-TAU_M-1);
    end
    
    % run KF from KFstart up until time of recent measurement
    for t_KF = KF_start:(t-TAU_M)
        
        if(cov_prior_adj)
                        
            % determine tNoACK vector for specific filter step
            if(N_VEH==1)
                
                if( (t_KF+TAU_AC-1)>0 )
                    t_NoACK_KF = t_NoACK(1,t_KF+TAU_AC-1);
                else
                    t_NoACK_KF = 1;
                end
                
            else
                
                % one-step approx
                t_NoACK_KF = 1;
                
            end
            
            % prepare control options for use in cov. prior adj.
            % u_options is indexed backwards
            % u_options(:,1) is most recent, (:,tNoACK_KF) is furthest back
            u_options = zeros(N_CONTROLS_ALL,t_NoACK_KF);
            for k = 1:t_NoACK_KF
                if(t_KF-k<1)
                    bTMP = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,1);
                else
                    bTMP = b_no_loss(:,t_KF-k);
                end
                u_options(:,k) = E*M^k*bTMP;
            end
            
            % overwrite u_options to zero if channel is ACK'd
            if(N_VEH>1)
                for i = 1:N_VEH
                    if( (t_KF+TAU_AC(i)-1)>0 && ...
                            (t_KF+TAU_AC(i)-1)<SIM_LEN && ...
                            (t_NoACK(i,t_KF+TAU_AC(i)-1)==0) )
                        u_options(i,1) = 0;
                    end
                end
            end
            
            
        else
            
            u_options = [];
            t_NoACK_KF = [];
            
        end
        
        if(t_KF<=1)
            A_KF = eye(size(A_SYS));
            Dc_KF_hat = makeDc(zeros(N_VEH,1),zeros(N_VEH,1),N_CONTROLS_VEH,N_HORIZON);
            Xh_in = [X_HAT1;zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,1)];
            P_in = P_1;
            U_in = zeros(N_CONTROLS_VEH*N_HORIZON*N_VEH,1);
            y_in = yh(:,1);
            Dm_in = Dm(:,:,1);
        else
            A_KF = A_SYS;
            Dc_KF_hat = Dc_hat(:,:,t_KF-1);
            Xh_in = Xh(:,t_KF-1);
            P_in = P(:,:,t_KF-1);
            U_in = U(:,t_KF-1);
            y_in = yh(:,t_KF);
            Dm_in = Dm(:,:,t_KF);
        end
        
        if(cov_prior_adj)
            cov.u_options = u_options;
            cov.t_NoACK = t_NoACK_KF;
        else
            cov=0;
        end
        
        if(print_debug_KF~=0)
            pd_KF.t = t; pd_KF.t_KF = t_KF;
        else
            pd_KF = 0;
        end
            
        % Xh(:,t-tm): xHat_{t-tm|t-tm},bHat_{t-tm-1}
        [Xh(:,t_KF),P(:,:,t_KF),Pstar_save_out] = JLSKF(Xh_in,P_in,y_in,U_in,Dc_KF_hat,...
            NX_SYS,N_VEH,N_CONTROLS_VEH,N_HORIZON,Dm_in,A_KF,Bu_SYS,E,M,C_SYS,...
            W_KF,V_KF,ALPHAC_BAR,cov,pd_KF);
        
        % log Pstar 
        if(logPstar)
            
            if(isempty(Pstar_save{t_KF}))
                
                % log Pstar array 
                Pstar_save{t_KF} = Pstar_save_out;
                
            else
                
                % backup-and-rerun is overwriting previously set Pstar
                if(print_debug_KF)
                    fprintf('\nt=%d, KF tKF=%d, overwriting Pstar\n',t,t_KF)
                end
                
                if(N_VEH==1)
                    % If delays are not excessively long vs. n_ACKHistory,
                    %   then only overwrites should be with empty Pstar
                    % --> keep the original Pstars in Pstar_save, 
                    %   and indicate that this step of Pstar is overwritted
                    if(isempty(Pstar_save_out))
                        Pstar_overwrite(t_KF) = 1;
                    else
                        disp('warning -- OVERWRITE WITH NONEMPTY PSTAR')
                        disp(Pstar_save_out)
                    end
                else
                    % makes a cell with each successive Pstar overwrite.  
                    tmp = Pstar_save{t_KF};
                    if(iscell(tmp))
                        tmpsize = length(tmp);
                        tmpcell = tmp;
                    else
                        tmpsize = 1;
                        tmpcell{1} = tmp;
                    end
                    PstarHistory = cell(1,tmpsize+1);
                    PstarHistory(1:tmpsize) = tmpcell;
                    PstarHistory{tmpsize+1} = Pstar_save_out;
                    Pstar_save{t_KF} = PstarHistory;
                    clear tmpcell
                end
            end
        end
        
        if(print_debug)
            if(cov_prior_adj)
                if(N_VEH==1)
                    fprintf('\nt=%d, KF tKF=%d, tNoACK_KF(%d)=%d\n',...
                        t,t_KF,t_KF+TAU_AC-1,t_NoACK_KF)
                else
                    %   fprintf('\nt=%d, KF td=%d, tNoACK_KF(%d)=%d\n',...
                    %      t,td,td+tac(i)-1,tNoACK_KF(i))
                    % (add in once resolve ambiguity)
                end
            end
            if(size(A_SYS,1)==1)
                % only print estimate for scalar sys
                disp_start = t_KF - 5;
                if(disp_start<1); disp_start = 1; end
                fprintf('Xh(:,%d:%d)=\n',disp_start,t_KF);...
                    disp(Xh(:,disp_start:t_KF)')
            end
        end
        
    end
    
    %%%%%%%
    % if control is to be computed/sent this step
    if(max(PI_C(:,t)))
        
        %%%%%%%
        
        MPCtic = tic;
        
        % grab constraints (to accommodate time-varying)
        [u_max,u_min,xmax,xmin] = paramsNow(U_MAX_T,U_MIN_T,X_MAX,X_MIN,...
            t+1,N_HORIZON);
        
        % Forward propagation: XhMPC and k_p^i's
        % starts with Xh: xHat_{t-tm|t-tm}, bHat_{t-tm-1}
        % compute xHatMPC_{t-tm+1:t+tc|t-tm}, bHatMPC_{t-tm:t+tc-1}
        % first step: uses uHat_{t-tm} <-- Dh(:,:,t-tm)
        % *DEBUG NOTE* Goal of XhMPC(:,end,:) is to match true X(:,:)
        
        U_fwd = U(:,(t-TAU_M):(t+TAU_C-1));
        D_fwd = Dc_hat(:,:,(t-TAU_M):(t+TAU_C-1));
        [XhMPC(:,:,t+TAU_C),p_i] = prepMPC(t,Xh(:,t-TAU_M),U_fwd,D_fwd,...
            PI_C,A_SYS,Bu_SYS,E,M,NX_SYS,N_VEH,N_CONTROLS_VEH,N_HORIZON,TAU_M,TAU_C);
              
        solve_status = 0;
        counter = 1;
        T_MPC = N_HORIZON;
        while(solve_status==0)
            
            % note -- logged variables corresponding to U(t+TAU_C) have
            % time index t+TAU_C as well 
            % (even though control computed at time t)
            
            % compute U_{t+tc}^i, forall i s.t. {Pi_c(i,t) = 1}
            [U_MPC,Jcomp(t+TAU_C),status,X_plan_out,~] = schedMPC(XhMPC(1:NX_SYS,...
                end,t+TAU_C),XhMPC((NX_SYS+1):end,end,t+TAU_C),p_i,...
                T_MPC,A_SYS,Bu_SYS,M,E,QMPC,QfMPC,RMPC,u_max,u_min,...
                xmin,xmax,[]);
            
            if(strfind(status,'Solved'))
                solve_status=1;
                if(print_debug)
                    fprintf('\nt=%d, MPC: %s\n',t,status)
                end
                
            elseif( strcmp(status,'Failed') )
                disp('FAILED')
                
            elseif( strcmp(status,'Infeasible') )
                disp('INFEASIBLE')
                disp(counter)
                disp(XhMPC(:,end,t+TAU_C))
                T_MPC = T_MPC+4;
                [U_MAX,U_MIN,xmax,xmin] = paramsNow(U_MAX_T,U_MIN_T,...
                    X_MAX,X_MIN,t+1,T_MPC);
            end
            
            counter = counter+1;
            if(counter>2)
                disp('MAXCOUNTER')
                U_MPC = zeros(N_CONTROLS_ALL,N_HORIZON);
                MPC_fail(t+TAU_C) = 1;
                break
            end
            
        end
        U_MPC = U_MPC(:,1:N_HORIZON);    % truncate if TMPC>Np
        X_plan(:,:,t+TAU_C) = X_plan_out(:,1:N_HORIZON);
        
        MPC_time(t) = toc(MPCtic);
        
    else % (for saving - make clear not set)
        
        U_MPC = zeros(N_CONTROLS_ALL,N_HORIZON);
        
    end
    
    % translate (NU x Np Umpc) into buffer shape Nu x Np x Nv
    U_vec = reshape(U_MPC',[N_CONTROLS_VEH*N_HORIZON*N_VEH,1]);
    if(isempty(strfind(CODEBOOK,'none')))
        [~,U(:,t+TAU_C)] = quantiz(U_vec,PARTITION,CODEBOOK);
    elseif(strfind(CODEBOOK,'none'))
        U(:,t+TAU_C) = U_vec;
    end
    
    % true system propagation
    % reshape into MJLS form (for step t), uses D_t(pi(t-tc),alpha_c(t-tc))
    if(t<=TAU_C)
        Dc(:,:,t) = makeDc(zeros(1,N_VEH),zeros(1,N_VEH),N_CONTROLS_VEH,N_HORIZON);
    else
        Dc(:,:,t) = makeDc(PI_C(:,t-TAU_C),alpha_c(:,t-TAU_C),...
            N_CONTROLS_VEH,N_HORIZON);
    end
    I = eye(size(Dc(:,:,1)));
    AA = [A_SYS,Bu_SYS*E*M*(I-Dc(:,:,t));zeros(N_VEH*N_HORIZON*N_CONTROLS_VEH,...
        NX_SYS),M*(I-Dc(:,:,t))];
    BU = [Bu_SYS*E*Dc(:,:,t);Dc(:,:,t)];
    
    % propagate system x_{t+1} = f(D_t,x_t,U_t,w_t)
    X(:,t+1) = AA*X(:,t) + BU*U(:,t) + BW*w_t(:,t);
    
    % actual applied u_t for saving
    u(:,t) = E*X(NX_SYS+1:end,t+1);
    
    % measurement z_{t+1}:
    y(:,t+1) = C_SYS*X(1:NX_SYS,t+1) + v_t(:,t+1);
    
    loop_time(t) = toc(looptic);
    
end

xF = X(1:NX_SYS,end);
% compute "actual" cost
jj = u(:,1)'*RMPC*u(:,1);  % control at step 1
for t = 2:SIM_LEN
    % states at step 2 (affected by u(1)), through Ns
    jj = jj + X(1:NX_SYS,t)'*QMPC*X(1:NX_SYS,t) + u(:,t)'*RMPC*u(:,t);
end
% final state (affected by final u)
Jsim = jj + xF'*QMPC*xF;


% populate results struct for output:

results.X = X;
results.u = u;
results.U = U;
 
results.u_no_loss = u_no_loss;
results.b_no_loss = b_no_loss;
results.t_NoACK = t_NoACK;

results.Xh = Xh;
results.P = P;

results.t_NoACK_save = t_NoACK_save;
results.Pstar_save = Pstar_save;
if(N_VEH==1);results.Pstar_overwrite = Pstar_overwrite;end

results.XhMPC = XhMPC;
results.Jcomp = Jcomp;
results.X_plan = X_plan;
results.MPC_time = MPC_time;
results.MPC_fail = MPC_fail;

results.loop_time = loop_time;
results.Jsim = Jsim/SIM_LEN;
results.rms_est_error = nanrms(X(1:NX_SYS,:) - Xh(1:NX_SYS,:),2);
results.rms_state_error = nanrms(X(1:NX_SYS,:),2);

% run-specific parameters
results.P1 = P_1;
results.v_t = v_t;
results.w_t = w_t;
results.x_IC = X_IC;
results.X_HAT1 = X_HAT1;
results.alpha_c = alpha_c;
results.alpha_m = alpha_m;
results.alpha_a = alpha_a;
results.PI_C = PI_C;
results.PI_M = PI_M;
results.PI_A = PI_A;

% system
sys.A = A_SYS;
sys.Bu = Bu_SYS;
sys.C = C_SYS;
sys.U_MAX = U_MAX;
sys.U_MIN = U_MIN;
sys.Q = QMPC;
sys.Qf = QfMPC;
sys.R = RMPC;
sys.N_p = N_HORIZON;
sys.W = W_KF;
sys.V = V_KF;
sys.alphac_bar = ALPHAC_BAR;
sys.TAU_A = TAU_A;
sys.TAU_C = TAU_C;
sys.TAU_M = TAU_M;
sys.n_ACKHistory = N_ACKHISTORY;

results.sys = sys;

end





