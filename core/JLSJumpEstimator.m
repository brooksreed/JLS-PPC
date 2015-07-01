function [Dc_hat_mod,alphac_hat_mod,Da_mod,KF_start_out,...
    t_NoACK_mod,got_ACK] = JLSJumpEstimator(Dc_hat_mod,PI_C,...
    Da_mod,alpha_c,alphac_hat_mod,PI_A,T_S,alpha_a,t_sim,TAU_M,TAU_C,...
    TAU_A,TAU_AC,N_controls,N_p,t_NoACK_mod,N_ACKHISTORY,print_debug_JE)
% Jump estimator and prep for (lossy,delayed) KF
% [Dc_hat_mod,alphac_hat_mod,Da_mod,KF_start_out,...
%     t_NoACK_mod,got_ACK] = JLSJumpEstimator(Dc_hat_mod,PI_C,...
%     Da_mod,alpha_c,alphac_hat_mod,PI_A,T_S,alpha_a,t_sim,TAU_M,TAU_C,...
%     TAU_A,TAU_AC,N_controls,N_p,t_NoACK_mod,N_ACKHISTORY,print_debug_JE)
% Updates Dc_hat, alphac_hat based on delayed ACKs
% Determines time that ACK'd buffer starts, outputs KF_start
% Keeps track of t_NoACK counter, also outputs gotACK flag
%
% Algorithm will back up as far as it can towards the most recent time of
% an ACK, using ACK histories.  
% It does not back up extra far (does not use all ACK history if not
% needed).  
% t_NoACK(t-TAU_A) resets to zero if ACK received at time t 


n_ACKs = size(t_NoACK_mod,1);

% Algorithm uses previous step counter when figuring out how far to back up
if(t_sim-TAU_A-1>0)
    t_NoACK_JE = t_NoACK_mod(:,t_sim-TAU_A-1);
else
    t_NoACK_JE = zeros(n_ACKs);
end

% limit "lookback" to the length of the ACK history sent
% "-1" is so that no *extra* backups done if 1 ACK sent (nACKHistory=1)
for i = 1:n_ACKs
    if(t_NoACK_JE(i) >= N_ACKHISTORY)
        t_NoACK_JE(i) = N_ACKHISTORY-1;
    end
end

% determine ACKs available at this step
if(t_sim-TAU_A<1)
    a_inds = [];
else
    Da_mod(:,:,t_sim-TAU_A) = diag(...
        PI_A(:,t_sim-TAU_A).*alpha_a(:,t_sim-TAU_A));
    a_inds = find(diag(Da_mod(:,:,t_sim-TAU_A)));
end

% use ACK information to update Dhat and alphaHat

% check at start - don't reference negative times
% slight hack, more precise would check individual channels... 
check_start = (t_sim-TAU_A-max(TAU_AC)>TAU_C);
if( ~isempty(a_inds) && check_start )
    
    % for each ACK channel
    max_t_ACK = zeros(size(a_inds));min_t_ACK = NaN*zeros(size(a_inds));
    t_ACK = cell(1,length(a_inds));
    for i = a_inds
        
        % find time range for new ACKs
        % (furthest back in time ACK history will be useful)
        t_ACK{i} = t_sim-TAU_A-TAU_AC(i) - ...
            t_NoACK_JE(i):t_sim-TAU_A-TAU_AC(i);
        
        % constrain tNoACKAlg to be nACKHistory-1?
        % or increment left side +1
        
        t_ACK{i} = t_ACK{i}(t_ACK{i}>TAU_C);
        min_t_ACK(i) = min(t_ACK{i});
        max_t_ACK(i) = max(t_ACK{i});
    end
    
    % back-up each ACK channel to appropriate time
    for i = a_inds
        % run fwd incorporating new ACKs
        for tprime = t_ACK{i}
            % update specific alpha_cHat's
            alphac_hat_mod(i,tprime-TAU_C) = ...
                alpha_c(i,tprime-TAU_C);
        end
    end
    
    % backup to earliest time across all channels
    backup_start = min(min_t_ACK);
    backup_end = max(max_t_ACK);
    
    % update Dhat using alphahat
    for tprime = backup_start:backup_end
        Dc_hat_mod(:,:,tprime) = makeD_c(PI_C(:,tprime-TAU_C),...
            alphac_hat_mod(:,tprime-TAU_C),N_controls,N_p);
    end
    
    % KFstart is the oldest (a posteriori) step of backup-rerun
    % KFstart uses the oldest updated Dhat(backupStart)
    KF_start_out = min((t_sim-TAU_M),(backup_start+1));
    if(print_debug_JE)
       fprintf('\nt=%d, JE t-tm=%d, backupStart+1=%d\n',...
           t_sim,t_sim-TAU_M,backup_start+1) 
    end
    
else
    
    KF_start_out = t_sim - TAU_M;   % usual 1-step update
    
end

% update gotACK flag for each channel for KF to use
got_ACK = zeros(size(Da_mod,3));
for i = 1:length(a_inds)
    got_ACK(a_inds(i)) = 1;
    if(print_debug_JE)
        fprintf('\nt=%d, JE GOT ACK, channel %d \n',t_sim,i)
    end
end

% increment tNoACK (plus lookahead)
if( (t_sim-TAU_A-1)>0 )
    
    nL = 100;
    
    % some checks for the end of the mission 
    if(t_sim-TAU_A+nL>size(t_NoACK_mod,2))
        max_add = size(t_NoACK_mod,2);
    else
        max_add = t_sim-TAU_A+nL;
    end
    lookahead = 1:(max_add-(t_sim-TAU_A));
    
    % update tNoACK based on new ACKs this step (t-ta)
    for i = 1:n_ACKs
        if(got_ACK(i))
            t_NoACK_mod(i,t_sim-TAU_A) = 0;
            
            % determine history able to be updated based on this ACK
            % () for clarity - (if nACKHistory=1: just t-ta)
            t_backup = (t_sim-TAU_A)-(N_ACKHISTORY-1);   
            if(t_backup<1)
                t_backup = 1;
            end
            
            % determine future that is known based on this ACK
            % (due to schedule -- time until next planned ACK)
            
            % (note - this assumes one control/ack pair per schedule)
            %   tau_ac implicitly used for getting ACK
            %       fwd knowledge based entirely on Ts 
            %       Ts into future, new ACK is expected (counter starts up)
            t_fwd = (t_sim-TAU_A)+(T_S-1);
            
            % zero counter based on ACKs received
            t_NoACK_mod(i,t_backup:t_fwd) = 0; 
            
            % shift lookahead according to schedule
            t_NoACK_mod(i,t_fwd+1:max_add) = lookahead(1:(max_add-t_fwd));
            
        end
    end
end

if(print_debug_JE && (t_sim-TAU_A>1) )
    for i = 1:n_ACKs
        fprintf('\nt=%d, JE tNoACK(t-ta)=%d, tNoACK^%d(t)=%d\n',...
            t_sim,t_NoACK_mod(i,t_sim-TAU_A),...
            i,t_NoACK_mod(i,t_sim))
    end
    if(exist('backupStart','var'))
        fprintf('      JE backupStart: t=%d\n',backup_start)
    end
end

if(print_debug_JE)
    for i = 1:n_ACKs
        if(t_NoACK_mod(i) > N_ACKHISTORY)
            fprintf('\nt=%d, JE WARNING: #dropped ACKs > ACK history, channel %d \n',t_sim,i)
        end
    end
end


