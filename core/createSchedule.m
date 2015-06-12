function [Pi_c,xi,lambda,tap,Ts] = createSchedule(sched,Nv,Ns,tc)
% creates Pi_c, xi, lambda schedule
% [Pi_c,xi,lambda,tap,Ts] = createSchedule(sched,Nv,Ns,tc)
% Nv: #vehicles, Ns: sim length, tc: control delay
% sched options:
% 'RRmeas' -- round-robin for meas, with 1 unused slot (match Ts of MX)
% 'MX_noACK', 'IL_noACK' 
% 'MX_piggyback', 'IL_piggyback' -- ACK sched matches measurements
% 'MX_ACK1', 'IL_ACK1' -- ACK 1 slot after control (assumes tc=1)
%
% 'SISO_' options for _piggyback or _noACK:
% 'SISOALL' - [1],[1] (std. discrete-time)
% 'SISO2' - [1 0], [0 1] for Pi_c, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for Pi_c, xi
%
% tap = Nv x 1 vector of time between planned control RX and ACK TX

% BR, 4/23/2014

% update 6/17/2014 for ACKs
% kind of sloppy right now, lots of hardcoding...  make more modular?

if(strfind(sched,'SISOALL'))
    Ts = 1;
    Pi_cBase = 1;
    xiBase = 1;
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = 0;
    end
end

if(strfind(sched,'SISO2'))
    Ts = 2;
    Pi_cBase = [1 0];
    xiBase = [0 1];
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = [0 0];
    end
end

if(strfind(sched,'SISO2ALLCONTROL'))
    Ts = 2;
    Pi_cBase = [1 1];
    xiBase = [0 1];
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = [0 0];
    end
end

if(strfind(sched,'SISO4'))
    Ts = 4;
    Pi_cBase = [1 0 0 0];
    xiBase = [0 0 1 0];
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = [0 0 0 0];
    end
end


if(strfind(sched,'RRmeas'))
    sched = 'RRmeas_noACK';
    Ts = Nv+1;  % set RR sched equal to MX sched length
    Pi_cBase = ones(Nv,Ts);
    xiBase = zeros(Nv,Ts);
    for i = 1:Nv
        xiBase(i,i) = 1;
    end
    lambdaBase = zeros(Nv,Ts);
end

if( ~isempty((strfind(sched,'noACK'))) ||...
        ~isempty((strfind(sched,'piggyback'))) )

    if(strfind(sched,'block'))

        Ts = Nv + Nv;
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,Nv+i) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,i) = 1;
        end
        
    elseif(strfind(sched,'IL'))

        Ts = Nv + Nv;
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,(2*i)) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,(2*i-1)) = 1;
        end
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + 1;
        Pi_cBase = zeros(Nv,Ts);
        Pi_cBase(:,Ts) = ones(Nv,1);
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,i) = 1;
        end
        
    end
    
    % set lambda
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    else
        lambdaBase = zeros(size(xiBase));
    end

elseif(strfind(sched,'ACK1'))   % dedicated ACK slot right after u
    
    if(strfind(sched,'IL'))
        Ts = Nv + Nv + Nv;
        lambdaBase = zeros(Nv,Ts);
        for i = 1:Nv
            lambdaBase(i,(3*i)) = 1;
        end        
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,(3*i-1)) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,(3*i-2)) = 1;
        end        
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + Nv + 1;
        
        Pi_cBase = zeros(Nv,Ts);
        Pi_cBase(:,Nv+1) = ones(Nv,1);
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,i) = 1;
        end        
        lambdaBase = zeros(Nv,Ts);
        for i = 1:Nv
            lambdaBase(i,Nv+1+i) = 1;
        end
        
    end
    
end



% compute tap
tap = zeros(Nv,1);
if(isempty(strfind(sched,'noACK')))
    for i = 1:Nv
        Pi_cPos = find(Pi_cBase(i,:))+tc;   % planned control RX
        lambdaPos = find(lambdaBase(i,:));
        if(Pi_cPos>lambdaPos)
            lambdaPos = lambdaPos + Ts;
        end
        try
            tap(i) = lambdaPos - Pi_cPos;
        catch
            % (SISO2ALLCONTROL -- fix eventually)
            tap(i) = 0;
        end
    end
end

nPers = ceil(Ns/Ts);
Pi_c = repmat(Pi_cBase,1,nPers);
xi = repmat(xiBase,1,nPers);
lambda = repmat(lambdaBase,1,nPers);
Pi_c = Pi_c(:,1:Ns);
xi = xi(:,1:Ns);
lambda = lambda(:,1:Ns);

