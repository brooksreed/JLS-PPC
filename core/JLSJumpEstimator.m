function [D_cHat,alpha_cHat,D_a,KFstart,tNoACK] = JLSJumpEstimator(...
    D_cHat,Pi_c,D_a,alpha_c,alpha_cHat,Pi_a,alpha_a,t,tm,tc,ta,tac,...
    Nu,Np,tNoACK,nACKHistory,printDebug)
% Jump estimator and prep for (lossy,delayed) KF
% [D_cHat,alpha_cHat,D_a,KFstart,tNoACK] = JLSJumpEstimator(...
%    D_cHat,Pi_c,D_a,alpha_c,alpha_cHat,Pi_a,alpha_a,t,tm,tc,ta,tac,...
%    Nu,Np,tNoACK,nACKHistory,printDebug)
% Updates D_cHat based on delayed ACKs
% Determines time that ACK'd buffer starts, outputs KFstart
%


% BR, 5/27/2015 (modified from single ACK version)
% 5/28 - ACK histories


% note -- actual sending of ack histories can be compressed a lot if
% consider control schedule (pi), just needs to be re-expanded on the
% estimator end
% - future: think about how would work in exps

% initial approach - clarity over efficiency
% refactor later?
% update help @ top...


% Algorithm will back up as far as it can towards the most recent time of
% an ACK, using ACK histories.  
% It does not back up extra far (does not use all ACK history if not
% needed).  

% tNoACK counter is updated with a "lag" of ta
% tNoACK(t-ta) resets to zero if ACK received at time t 

nACKs = size(tNoACK,1);

% Algorithm uses previous step counter when figuring out how far to back up
% tNoACK for step t-ta-1 used by algorithm
if(t-ta-1>0)
    tNoACKAlg = tNoACK(:,t-ta-1);
else
    tNoACKAlg = zeros(nACKs);
end
% limit "lookback" to the length of the ACK history sent
for i = 1:nACKs
    if(tNoACKAlg(i) > nACKHistory)
        tNoACKAlg(i) = nACKHistory;
    end
end

% determine ACKs available at this step
% at step t, ACKs are sent at t-ta
if(t-ta<1)
    aInds = [];
else
    D_a(:,:,t-ta) = diag(Pi_a(:,t-ta).*alpha_a(:,t-ta));
    aInds = find(diag(D_a(:,:,t-ta)));
end

% use ACK information to update Dhat and alphaHat

% check at start - don't reference negative times
% slight hack, more precise check would individual channels... 
checkStart = (t-ta-max(tac)>tc);
if( ~isempty(aInds) && checkStart )
    
    % for each ACK channel
    min_tACK = zeros(size(aInds));max_tACK=min_tACK;
    tACK = cell(1,length(aInds));
    for i = aInds
        % find time range for new ACKs
        % (furthest back in time ACK history will be useful)
        tACK{i} = t-ta-tac(i)-tNoACKAlg(i):t-ta-tac(i);
        tACK{i} = tACK{i}(tACK{i}>tc);
        min_tACK(i) = min(tACK{i});
        max_tACK(i) = max(tACK{i});
    end
    
    % back-up each ACK channel to appropriate time
    for i = aInds
        % run fwd incorporating new ACKs
        for tprime = tACK{i}
            % update specific alphaHat's
            alpha_cHat(i,tprime-tc) = alpha_c(i,tprime-tc);
        end
    end
    
    % backup to earliest time across all channels
    backupStart = min(min_tACK);
    backupEnd = max(max_tACK);
    
    % update Dhat using alphahat
    for tprime = backupStart:backupEnd
        D_cHat(:,:,tprime) = makeD_c(Pi_c(:,tprime-tc),alpha_cHat(:,tprime-tc),...
            Nu,Np);
    end
    
    % KFstart is the oldest (a posteriori) step of backup-rerun
    % KFstart uses the oldest updated Dhat(backupStart)
    KFstart = min((t-tm),(backupStart+1));
    if(printDebug)
       fprintf('\nt=%d, JE t-tm=%d, backupStart+1=%d\n',t,t-tm,backupStart+1) 
    end
    
else
    
    KFstart = t - tm;   % usual 1-step update
    
end

% update gotACK flag for each channel for KF to use
gotACK = zeros(size(D_a,3));
for i = 1:length(aInds)
    gotACK(aInds(i)) = 1;
    if(printDebug)
        fprintf('\nt=%d, JE GOT ACK, channel %d \n',t,i)
    end
end

% increment tNoACK (plus lookahead)
if( (t-ta-1)>0 )
    % (once this works, modify for multivar)
    
    nL = 10;
    
    % some checks for the end of the mission 
    if(t-ta+nL>length(tNoACK))
        maxadd = length(tNoACK);
    else
        maxadd = t-ta+nL;
    end
    lookahead = 1:(maxadd-(t-ta));
    
    % increment future starting from current value
    %startVal = tNoACK(:,t-ta-1);
    %tNoACK(:,(t-ta):maxadd)  =  lookahead + startVal;
    
    % update tNoACK based on new ACKs this step (t-ta)
    for i = 1:nACKs
        if(gotACK(i))
            tNoACK(i,t-ta) = 0;
            
            % zero past counters based on history
            tBackup = t-nACKHistory-ta;
            if(tBackup<1)
                tBackup = 1;
            end
            tNoACK(i,tBackup:t-ta) = 0;
            
            % increment future, starting at 1
            tNoACK(i,(t-ta+1):maxadd) = lookahead;
                        
        else
            
            startVal = tNoACK(i,t-ta-1);
            tNoACK(i,(t-ta):maxadd-1) = startVal + lookahead;
            
        end
    end
end

if(printDebug && (t-ta>1) )
    for i = 1:nACKs
        fprintf('\nt=%d, JE tNoACK(t-ta)=%d, tNoACK(t)=%d\n',t,tNoACK(i,t-ta),tNoACK(i,t))
    end
    if(exist('backupStart','var'))
        fprintf('    JE furthest back posterior estimate: t=%d\n',backupStart)
    end
end

if(printDebug)
    for i = 1:nACKs
        if(tNoACK(i) > nACKHistory)
            fprintf('\nt=%d, JE WARNING: #dropped ACKs > ACK history, channel %d \n',t,i)
        end
    end
end


