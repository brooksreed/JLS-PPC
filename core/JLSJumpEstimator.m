function [Dh,alphat,a,KFstart,tNoACK] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np,tNoACK,nACKHistory)
% Jump estimator and prep for (lossy,delayed) KF
% [Dh,alphat,a,KFstart] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
%    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np)
% Updates Dh based on delayed ACKs
% (If no ACK, Dh is left as input -- currently uses alphaBar)
% Determines time that ACK'd buffer starts, outputs KFstart

% BR, 5/27/2015 (modified from single ACK version)
% 5/28 - ACK histories

% note -- actual sending of ack histories can be compressed a lot if
% consider control schedule (pi), just needs to be re-expanded on the
% estimator end
% - future: think about how would work in exps

% initial approach - clarity over efficiency
% refactor later?

printDebug = 1;

tNoACK  = tNoACK+1; % always increment 1 step at beginning
if(printDebug)
    fprintf('\ntNoACK=%d\n',tNoACK)
end

% limit "lookback" to the length of the ACK history sent
for i = 1:length(tNoACK)
    if(tNoACK(i) > nACKHistory)
        tNoACK(i) = nACKHistory;
        if(printDebug)
            fprintf('\n Step %d WARNING: #dropped ACKs > ACK history\n',t)
        end
    end
end

% determine ACKs available at this step
% at step t, ACKs are sent at t-ta
a(:,:,t-ta) = diag(Lambda(:,t-ta).*gamma(:,t-ta));
aInds = find(diag(a(:,:,t-ta)));

% check at startup - ack referencing negative time
if(~isempty(aInds))
    initTimes = tNoACK(aInds);
end

% run through algorithm if useful ACK has arrived
if( ~isempty(aInds) && (t-ta-(max(initTimes)))>1 )
    
    % for each ACK channel
    min_tACK = zeros(size(aInds));max_tACK=min_tACK;
    for i = aInds
        % find time range for new ACKs
        tACK{i} = t-ta-tap(i)-tNoACK(i):t-ta-tap(i);
        tACK{i} = tACK{i}(tACK{i}>tc);
        min_tACK(i) = min(tACK{i});
        max_tACK(i) = max(tACK{i});
    end
    
    % back-up each ACK channel to appropriate time
    for i = aInds
        % run fwd incorporating new ACKs
        for tprime = tACK{i}
            % update specific alphat's
            alphat(i,tprime-tc) = alpha(i,tprime-tc);
        end
    end
    
    % backup to earliest time across all channels
    backupStart = min(min_tACK);
    backupEnd = max(max_tACK);
    
    % update Dhat using alphahat
    for tprime = backupStart:backupEnd
        Dh(:,:,tprime) = makeD(Pi(:,tprime-tc),alphat(:,tprime-tc),...
            Nu,Np);
    end
    
    % if ACKs, KF start goes back to oldest useful (eg new) ACK
    KFstart = min((t-tm),backupStart);
    
else
    
    KFstart = t - tm;   % usual 1-step update
    
end

% reset ACK counter to zero for ACK channels received
for i = 1:length(aInds)
    tNoACK(aInds(i)) = 0; 
    if(printDebug)
        fprintf('\n Step %d GOT ACK\n',t)
    end
end
