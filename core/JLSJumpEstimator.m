function [Dh,alphaHat,a,KFstart,gotACK,tNoACK] = JLSJumpEstimator(Dh,Pi,a,alpha,alphaHat,...
    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np,tNoACK,nACKHistory,printDebug)
% Jump estimator and prep for (lossy,delayed) KF
% [Dh,alphaHat,a,KFstart] = JLSJumpEstimator(Dh,Pi,a,alpha,alphaHat,...
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

% Algorithm will back up as far as it can towards the most recent time of
% an ACK, using ACK histories.  
% It does not back up extra far (does not use all ACK history if not
% needed).  

% tNoACK counter is updated with a "lag" of ta
% tNoACK(t-ta) resets to zero if ACK received at time t 

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
    a(:,:,t-ta) = diag(Lambda(:,t-ta).*gamma(:,t-ta));
    aInds = find(diag(a(:,:,t-ta)));
end


% if useful ACK has arrived, update Dhat
% using alphaHat (alpha as known at estimator)

% check at start - don't reference negative times
% slight hack, more precise check would individual channels... 
checkStart = (t-ta-max(tap)>tc);
if( ~isempty(aInds) && checkStart )
    
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
            % update specific alphaHat's
            alphaHat(i,tprime-tc) = alpha(i,tprime-tc);
        end
    end
    
    % backup to earliest time across all channels
    backupStart = min(min_tACK);
    backupEnd = max(max_tACK);
    
    % update Dhat using alphahat
    for tprime = backupStart:backupEnd
        Dh(:,:,tprime) = makeD(Pi(:,tprime-tc),alphaHat(:,tprime-tc),...
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
gotACK = zeros(size(a,3));
for i = 1:length(aInds)
    gotACK(aInds(i)) = 1;
    if(printDebug)
        fprintf('\nt=%d, JE GOT ACK, channel %d \n',t,i)
        % only going to backup if got ACK...
    end
end

% increment true tNoACK (into the future)
if( (t-ta-1)>0 )
    
    % depending on delays, schedule of ACKs vs. measurements, the KF may
    % run further ahead than the ACKs have been updated for
    % rather than implement a separate counter inside KF, we fill in a
    % limited lookahead of tNoACK increments 
    lookahead = 10;
    
    % (once this works, modify for multivar)
    
    % some checks for the end of the mission 
    if(t-ta+lookahead>length(tNoACK))
        maxadd = length(tNoACK);
    else
        maxadd = t-ta+lookahead;
    end
    startVal = tNoACK(:,t-ta-1);
    %newVec = ones(1,maxadd - (t-ta)+1);
    %newVec1 = 0:(maxadd - (t-ta));
    newVec1 = 1:(maxadd-(t-ta)+1);
    % increment future starting from current value
    tNoACK(:,(t-ta):maxadd)  =  newVec1 + startVal;
    
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
            newVec2 = 1:(maxadd - (t-ta));
            tNoACK(i,(t-ta+1):maxadd) = newVec2;
                        
        end
    end
end

if(printDebug && (t-ta>1) )
    for i = 1:nACKs
        fprintf('\nt=%d, JE tNoACK(t-ta)=%d, tNoACK(t)=%d\n',t,tNoACK(t-ta),tNoACK(i,t))
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


