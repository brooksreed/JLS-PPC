function [pi,xi,lambda,tap,Ts] = createSchedule(sched,Nv,Ns,tc)
% creates pi, xi, lambda schedule
% [pi,xi,lambda,tap,Ts] = createSchedule(sched,Nv,Ns,tc)
% Nv: #vehicles, Ns: sim length, tc: control delay
% sched options:
% 'RRmeas' -- round-robin for meas, with 1 unused slot (match Ts of MX)
% 'MX_noACK', 'IL_noACK' 
% 'MX_piggyback', 'IL_piggyback' -- ACK sched matches measurements
% 'MX_ACK1', 'IL_ACK1' -- ACK 1 slot after control (assumes tc=1)
%
% 'SISO_' options for _piggyback or _noACK:
% 'SISOALL' - [1],[1] (std. discrete-time)
% 'SISO2' - [1 0], [0 1] for pi, xi
% 'SISO4' - [1 0 0 0], [0 0 1 0] for pi, xi
%
% tap = Nv x 1 vector of time between planned control RX and ACK TX

% BR, 4/23/2014

% update 6/17/2014 for ACKs
% kind of sloppy right now, lots of hardcoding...  make more modular?

if(strfind(sched,'SISOALL'))
    Ts = 1;
    piBase = 1;
    xiBase = 1;
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = 0;
    end
end

if(strfind(sched,'SISO2'))
    Ts = 2;
    piBase = [1 0];
    xiBase = [0 1];
    if(strfind(sched,'piggyback'))
        lambdaBase = xiBase;
    elseif(strfind(sched,'noACK'))
        lambdaBase = [0 0];
    end
end


if(strfind(sched,'SISO4'))
    Ts = 4;
    piBase = [1 0 0 0];
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
    piBase = ones(Nv,Ts);
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
        piBase = zeros(Nv,Ts);
        for i = 1:Nv
            piBase(i,Nv+i) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,i) = 1;
        end
        
    elseif(strfind(sched,'IL'))

        Ts = Nv + Nv;
        piBase = zeros(Nv,Ts);
        for i = 1:Nv
            piBase(i,(2*i)) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,(2*i-1)) = 1;
        end
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + 1;
        piBase = zeros(Nv,Ts);
        piBase(:,Ts) = ones(Nv,1);
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
        piBase = zeros(Nv,Ts);
        for i = 1:Nv
            piBase(i,(3*i-1)) = 1;
        end
        xiBase = zeros(Nv,Ts);
        for i = 1:Nv
            xiBase(i,(3*i-2)) = 1;
        end        
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + Nv + 1;
        
        piBase = zeros(Nv,Ts);
        piBase(:,Nv+1) = ones(Nv,1);
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
        piPos = find(piBase(i,:))+tc;   % planned control RX
        lambdaPos = find(lambdaBase(i,:));
        if(piPos>lambdaPos)
            lambdaPos = lambdaPos + Ts;
        end
        tap(i) = lambdaPos - piPos;
    end
end

nPers = ceil(Ns/Ts);
pi = repmat(piBase,1,nPers);
xi = repmat(xiBase,1,nPers);
lambda = repmat(lambdaBase,1,nPers);
pi = pi(:,1:Ns);
xi = xi(:,1:Ns);
lambda = lambda(:,1:Ns);

