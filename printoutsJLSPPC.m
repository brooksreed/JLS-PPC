% run inside the loop of simMJLSPPC
% print out debug/track info

% BR, 5/20/2014

% add a "snapshot" plot feature?

% time (in code, at estimator) for printouts
tP = t;

fprintf('--------------------\nStep %d\n--------------------\n',tP)

fprintf('\nMeasurements received at %d:',tP)
disp(diag(S(:,:,tP-tm))')
fprintf('ACKs received at %d:',tP)
disp(diag(a(:,:,tP-ta))')

% (Need to save ackInds and ackTimes cell for this)
%{
if(~isempty(ackInds{t-ta}))
    for k = ackTimes{t-ta}
        fprintf('\nstep %d true Pi.*alpha:\n',k-tc)
        disp(Pi(:,k-tc).*alpha(:,k-tc))
        fprintf('step %d Pi.*alphaHat:\n',k-tc)
        disp(Pi(:,k-tc).*alphat(:,k-tc))
    end
end
%}


fprintf('---------\nKF innovations %d to %d\n\n',KFstart-1,tP-tm)
for td = KFstart:(tP-tm)
    bPrior = reshape(Xh((Nx+1):end,td),[Np,Nu*Nv])';


    fprintf('True vehicle buffers at %d: (',td-1)
    fprintf('Pi_%d:',td-1-tc)
    fprintf(' %d',Pi(:,td-1-tc))
    fprintf('; Pi.*alpha for %d: ',td-1-tc)
    fprintf('%d ',Pi(:,td-1-tc).*alpha(:,td-1-tc))
    fprintf('):\n')
    disp(reshape(X((Nx+1):end,td),[Np,Nu*Nv])')
    
    fprintf('Estimated buffers at %d:  (',td-1)
    fprintf('Pi.*alpha_tilde for %d: ',td-1-tc)
    fprintf('%0.2f  -- ',Pi(:,td-1-tc).*alphat(:,td-1-tc))
    ACKstr = cell(Nv,1);
    for i = 1:Nv
        if( (td-1-tc+tap(i)+ta < Ns ) && (a(i,i,td-1-tc+tap(i)+ta)==1) )
            ACKstr{i} = 'ACK';
        else
            ACKstr{i} = 'NoACK';
        end
        fprintf('U^%d_%d: %s ',i,td-1-tc,ACKstr{i})
    end
    fprintf('):\n')
    disp(bPrior)
    
    % x_{t-tm|t-tm}
    fprintf('a posteriori estimate: \nxHat_{%d|%d},  x_{%d}\n',...
    td,td,td)
    disp([Xh(1:Nx,td),X(1:Nx,td)])
    
end
    
% compare with predicted state 
% available Xh is t-tm
% available XhMPC is {t-tm+1:t+tc|t-tm}

tComp = tP-tm-tc-1;
if((tComp-tm>=1) && (max(Pi(:,tComp))))
    % COMPUTED at tComp = t-tm-tc-1
    % xHat_{tComp - tm|tComp-tm} --> xHat_{tComp+tc|tComp-tc}
    fprintf('\n--------------------\n')
    fprintf('MPC comp at t-tm-tc-1 = %d, applied t-tm-1=%d\n',...
        tComp,tComp+tc)
    % tComp-tm -> tComp + tc
    % [XHat_{tComp-tm|tComp-tm},XFwd_{tComp-tm+1 -> tComp+tc}]
    
    if(tComp-tm-tc-1>0)
        fprintf('Pi.*alpha_tilde for %d -> %d: \n',...
            tComp-tm-tc-1,tComp-tm)
        %fprintf('%0.2f     ',Pi(:,tComp-tm-tc-1:tComp-tm).*...
        %    alphat(:,tComp-tm-tc-1:tComp-tm))
        disp(Pi(:,tComp-tm-tc-1:tComp-tm).*alphat(:,tComp-tm-tc-1:tComp-tm))
        %fprintf('\n')    
        fprintf('uFwd for %d -> %d: \n',tComp-tm-1,tComp+tc-1)
        disp(E1*XhMPC((Nx+1:end),:,tComp+tc))
    end
    fprintf('xFwd for %d -> %d: \n',tComp-tm,tComp+tc)
    disp(XhMPC(1:Nx,:,tComp+tc))
    fprintf('-------\n')
    fprintf('uHat for %d -> %d: \n',tComp-tm-1,tComp+tc-1)
    disp(E1*Xh((Nx+1:end),tComp-tm:tComp+tc))
    fprintf('xHat for %d -> %d:\n',tComp-tm,tComp+tc)
    disp(Xh(1:Nx,tComp-tm:tComp+tc))
    fprintf('-------\n')
    fprintf('u for %d -> %d: \n',tComp-tm-1,tComp+tc-1)
    disp(E1*X((Nx+1:end),tComp-tm:tComp+tc))
    fprintf('x for %d -> %d:\n',tComp-tm,tComp+tc)
    disp(X(1:Nx,tComp-tm:tComp+tc))
    garbage=1;
end
fprintf('\n')
