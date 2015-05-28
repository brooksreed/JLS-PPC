function [Dh,alphat,a,KFstart] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np)
% Jump estimator and prep for (lossy,delayed) KF
% [Dh,alphat,a,KFstart] = JLSJumpEstimator(Dh,Pi,a,alpha,alphat,...
%    Lambda,gamma,t,tm,tc,ta,tap,Nu,Np)
% Updates Dh based on delayed ACKs
% (If no ACK, Dh is left as input -- currently uses alphaBar)
% Determines time that ACK'd buffer starts, outputs KFstart

% BR, 9/30/2014 (modified into separate function)

% determine ACKs available at this step
% at step t, ACKs are sent at t-ta
a(:,:,t-ta) = diag(Lambda(:,t-ta).*gamma(:,t-ta));
aInds = find(diag(a(:,:,t-ta)));
aTimes = t - ta - tap(aInds);
okInds = logical((aTimes-tc)>1);
ackInds = aInds(okInds);
ackTimes = aTimes(okInds);

if(~isempty(ackInds))
    % times where ACK will update Dh
    % ACK RD'd: go back and update DHat for ackInds at ackTimes
    % (Dh init to pi.*alphaBar -- only update times where ACK changes)
    for i = 1:length(ackInds)
        ii = ackInds(i);
        alphat(ii,ackTimes(i)-tc) = alpha(ii,...
            ackTimes(i)-tc);
    end
    for i = 1:length(ackInds)
        k = ackTimes((i));
        Dh(:,:,k) = makeD(Pi(:,k-tc),alphat(:,k-tc),Nu,Np);
    end
end

if(isempty(ackInds))
    KFstart = t - tm;   % usual 1-step update
else
    KFstart = min((t-tm),min(ackTimes));
end

return

%% future -- with ACK histories

%pass tNoACK in/out, input nACKHistory
%(think about how this might change with exp version, where don't have
%access to full history of acks in the sim)

tNoAck  = tNoAck+1; % always increment 1 step at beginning

% limit "lookback" to the length of the ACK history sent
for i = 1:Nv
    if(tNoAck(i) > nACKHistory)
        tNoAck(i) = nACKHistory;
        disp('WARNING: number of dropped ACKs exceeds ACK history')
    end
end
    
if(~isempty(ackInds))
    
    % times where ACK will update Dh
    % ACK RD'd: go back and update DHat for ackInds at ackTimes
    % (Dh init to pi.*alphaBar -- only update times where ACK changes)
    
    for i = 1:length(ackInds)
        ii = ackInds(i);
        % note... ackTimes indexed same as ackInds
        % tNoACK indexed over 1,...,Nv (so use ii)
        for tprime = ackTimes(i)-tNoACK(ii):ackTimes(i)
            alphat(ii,tprime-tc) = alpha(ii,...
                tprime-tc);
        end
    end
    
    % separate Dh loop because do whole makeD for all i at once
    for i = 1:length(ackInds)
        for tprime = ackTimes(i)-tNoACK(ii):ackTimes(i)
            Dh(:,:,tprime) = makeD(Pi(:,tprime-tc),alphat(:,tprime-tc),...
                Nu,Np);
        end
    end
    
end

if(isempty(ackInds))
    KFstart = t - tm;   % usual 1-step update
else
    % if ACKs, KF start goes back to oldest useful (eg new) ACK
    KFstart = min((t-tm),min(ackTimes-tNoAck(ackInds)));
end

for i = 1:length(ackInds)
    tNoACK(ackInds(i)) = 0; % for ackInds, reset counter to zero    
end

