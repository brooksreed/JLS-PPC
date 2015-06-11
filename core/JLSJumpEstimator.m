function [Dh,alphat,a,KFstart,gotACK,tNoACK] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np,tNoACK,nACKHistory,printDebug)
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
% update help @ top...

nACKs = size(tNoACK,1);

% tNoACK for step t-1 used by algorithm
% limit "lookback" to the length of the ACK history sent
if((t-1)>ta)
    tNoACKAlg = tNoACK(:,t-ta-1);
else
    tNoACKAlg = zeros(nACKs);
end
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
    a(:,:,t-ta) = diag(Lambda(:,t-ta).*gamma(:,t-ta));
    aInds = find(diag(a(:,:,t-ta)));
end

% check at startup - ack referencing negative time
% if(~isempty(aInds))
%     initTimes = tNoACKAlg(aInds);
% end

% run through algorithm if useful ACK has arrived
% && (t-ta-(max(initTimes)))>1
if( ~isempty(aInds) && (t-ta-max(tap)>tc))
    
    % for each ACK channel
    min_tACK = zeros(size(aInds));max_tACK=min_tACK;
    tACK = cell(1,length(aInds));
    for i = aInds
        % find time range for new ACKs
        % (furthest back in time ACK history will be useful)
        tACK{i} = t-ta-tap(i)-tNoACKAlg(i):t-ta-tap(i);
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
    % KF start is the a posteriori step
    KFstart = min((t-tm),(backupStart));
    
else
    
    KFstart = t - tm;   % usual 1-step update
    
end

gotACK = zeros(size(a,3));
% update gotACK flag for each channel for KF to use
for i = 1:length(aInds)
    gotACK(aInds(i)) = 1;
    if(printDebug)
        fprintf('\nt=%d GOT ACK, channel %d \n',t,i)
        % only going to backup if got ACK...
    end
end


% increment true tNoACK (into the future)
if( (t-ta-1)>0)
    
    if(t-ta+10>length(tNoACK))
        maxadd = length(tNoACK);
    else
        maxadd = t-ta+10;
    end
    startVal = tNoACK(:,t-ta-1);
    newVec = ones(1,maxadd - (t-ta)+1);
    newVec(1,:) = 0:(maxadd - (t-ta));
    tNoACK(:,(t-ta):maxadd)  =  newVec + startVal;
    %tNoACK(:,t-ta+1) = tNoACK(:,t-ta)+1;
    
    % update tNoACKSave properly
    for i = 1:nACKs
        if(gotACK(i))
            % if Jump Estimator says got ACK
            tNoACK(i,t-ta) = 0;
            tBackup = t-nACKHistory-ta;
            if(tBackup<1)
                tBackup = 1;
            end
            % update nACKHistory to past
            tNoACK(i,tBackup:t-ta) = 0;
            % increment future
            tNoACK(i,t-ta+1:maxadd) = 1:(maxadd-(t-ta));
        end
    end
end

if(printDebug)
    for i = 1:nACKs
        fprintf('\nt=%d, Jump Estimator tNoACK=%d\n',t,tNoACK(i,t))
    end
    if(exist('backupStart','var'))
        fprintf('furthest back posterior estimate: t=%d\n',backupStart)
    end
end

if(printDebug)
    for i = 1:nACKs
        if(tNoACK(i) > nACKHistory)
            fprintf('\nt=%d WARNING: #dropped ACKs > ACK history, channel %d \n',t,i)
        end
    end
end


