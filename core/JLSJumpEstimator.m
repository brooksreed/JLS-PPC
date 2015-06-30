function [D_cHat_global,alpha_cHat_global,D_a_global,KFstart,...
    tNoACK_global,gotACK] = JLSJumpEstimator(D_cHat_global,Pi_c_all,...
    D_a_global,alpha_c_all,alpha_cHat_global,Pi_a_all,Ts_sched,...
    alpha_a_all,t_global,tm_sched,tc_sched,ta_sched,tac_sched,...
    Ncontrols,Nhorizon,tNoACK_global,len_ACKHistory,printDebug_in)
% Jump estimator and prep for (lossy,delayed) KF
% [D_cHat_global,alpha_cHat_global,D_a_global,KFstart,...
%     tNoACK_global,gotACK] = JLSJumpEstimator(D_cHat_global,Pi_c_all,...
%     D_a_global,alpha_c_all,alpha_cHat_global,Pi_a_all,Ts_sched,...
%     alpha_a_all,t_global,tm_sched,tc_sched,ta_sched,tac_sched,...
%     Ncontrols,Nhorizon,tNoACK_global,len_ACKHistory,printDebug_in)
%
% Updates D_cHat, alpha_cHat based on delayed ACKs
% Determines time that ACK'd buffer starts, outputs KFstart
% Keeps track of tNoACK counter, also outputs gotACK flag
%
% Algorithm will back up as far as it can towards the most recent time of
% an ACK, using ACK histories.  
% It does not back up extra far (does not use all ACK history if not
% needed).  
% tNoACK(t-ta) resets to zero if ACK received at time t 


% BR, 5/27/2015 (modified from single ACK version)
% 5/28 - ACK histories
% v1.0 6/13/2015
% v1.1 6/16/2015

nACKs = size(tNoACK_global,1);

% Algorithm uses previous step counter when figuring out how far to back up
if(t_global-ta_sched-1>0)
    tNoACKAlg = tNoACK_global(:,t_global-ta_sched-1);
else
    tNoACKAlg = zeros(nACKs);
end

% limit "lookback" to the length of the ACK history sent
% "-1" is so that no *extra* backups done if 1 ACK sent (nACKHistory=1)
for i = 1:nACKs
    if(tNoACKAlg(i) >= len_ACKHistory)
        tNoACKAlg(i) = len_ACKHistory-1;
    end
end

% determine ACKs available at this step
if(t_global-ta_sched<1)
    aInds = [];
else
    D_a_global(:,:,t_global-ta_sched) = diag(...
        Pi_a_all(:,t_global-ta_sched).*alpha_a_all(:,t_global-ta_sched));
    aInds = find(diag(D_a_global(:,:,t_global-ta_sched)));
end

% use ACK information to update Dhat and alphaHat

% check at start - don't reference negative times
% slight hack, more precise would check individual channels... 
checkStart = (t_global-ta_sched-max(tac_sched)>tc_sched);
if( ~isempty(aInds) && checkStart )
    
    % for each ACK channel
    max_tACK = zeros(size(aInds));min_tACK = NaN*zeros(size(aInds));
    tACK = cell(1,length(aInds));
    for i = aInds
        
        % find time range for new ACKs
        % (furthest back in time ACK history will be useful)
        tACK{i} = t_global-ta_sched-tac_sched(i) - ...
            tNoACKAlg(i):t_global-ta_sched-tac_sched(i);
        
        % constrain tNoACKAlg to be nACKHistory-1?
        % or increment left side +1
        
        tACK{i} = tACK{i}(tACK{i}>tc_sched);
        min_tACK(i) = min(tACK{i});
        max_tACK(i) = max(tACK{i});
    end
    
    % back-up each ACK channel to appropriate time
    for i = aInds
        % run fwd incorporating new ACKs
        for tprime = tACK{i}
            % update specific alpha_cHat's
            alpha_cHat_global(i,tprime-tc_sched) = ...
                alpha_c_all(i,tprime-tc_sched);
        end
    end
    
    % backup to earliest time across all channels
    backupStart = min(min_tACK);
    backupEnd = max(max_tACK);
    
    % update Dhat using alphahat
    for tprime = backupStart:backupEnd
        D_cHat_global(:,:,tprime) = makeD_c(Pi_c_all(:,tprime-tc_sched),...
            alpha_cHat_global(:,tprime-tc_sched),Ncontrols,Nhorizon);
    end
    
    % KFstart is the oldest (a posteriori) step of backup-rerun
    % KFstart uses the oldest updated Dhat(backupStart)
    KFstart = min((t_global-tm_sched),(backupStart+1));
    if(printDebug_in)
       fprintf('\nt=%d, JE t-tm=%d, backupStart+1=%d\n',...
           t_global,t_global-tm_sched,backupStart+1) 
    end
    
else
    
    KFstart = t_global - tm_sched;   % usual 1-step update
    
end

% update gotACK flag for each channel for KF to use
gotACK = zeros(size(D_a_global,3));
for i = 1:length(aInds)
    gotACK(aInds(i)) = 1;
    if(printDebug_in)
        fprintf('\nt=%d, JE GOT ACK, channel %d \n',t_global,i)
    end
end

% increment tNoACK (plus lookahead)
if( (t_global-ta_sched-1)>0 )
    
    nL = 100;
    
    % some checks for the end of the mission 
    if(t_global-ta_sched+nL>size(tNoACK_global,2))
        maxadd = size(tNoACK_global,2);
    else
        maxadd = t_global-ta_sched+nL;
    end
    lookahead = 1:(maxadd-(t_global-ta_sched));
    
    % update tNoACK based on new ACKs this step (t-ta)
    for i = 1:nACKs
        if(gotACK(i))
            tNoACK_global(i,t_global-ta_sched) = 0;
            
            % determine history able to be updated based on this ACK
            % () for clarity - (if nACKHistory=1: just t-ta)
            tBackup = (t_global-ta_sched)-(len_ACKHistory-1);   
            if(tBackup<1)
                tBackup = 1;
            end
            
            % determine future that is known based on this ACK
            % (due to schedule -- time until next planned ACK)
            
            % (note - this assumes one control/ack pair per schedule)
            %   tau_ac implicitly used for getting ACK
            %       fwd knowledge based entirely on Ts 
            %       Ts into future, new ACK is expected (counter starts up)
            tFwd = (t_global-ta_sched)+(Ts_sched-1);
            
            % zero counter based on ACKs received
            tNoACK_global(i,tBackup:tFwd) = 0; 
            
            % shift lookahead according to schedule
            tNoACK_global(i,tFwd+1:maxadd) = lookahead(1:(maxadd-tFwd));
            
        end
    end
end

if(printDebug_in && (t_global-ta_sched>1) )
    for i = 1:nACKs
        fprintf('\nt=%d, JE tNoACK(t-ta)=%d, tNoACK^%d(t)=%d\n',...
            t_global,tNoACK_global(i,t_global-ta_sched),...
            i,tNoACK_global(i,t_global))
    end
    if(exist('backupStart','var'))
        fprintf('      JE backupStart: t=%d\n',backupStart)
    end
end

if(printDebug_in)
    for i = 1:nACKs
        if(tNoACK_global(i) > len_ACKHistory)
            fprintf('\nt=%d, JE WARNING: #dropped ACKs > ACK history, channel %d \n',t_global,i)
        end
    end
end


