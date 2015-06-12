function [Pi_c,Pi_m,Pi_a,tac,Ts] = createSchedule(sched,Nv,Ns,tc)
% creates Pi_c, Pi_m, Pi_a schedule
% [Pi_c,Pi_m,Pi_a,ta_c,Ts] = createSchedule(sched,Nv,Ns,tc)
% Nv: #vehicles, Ns: sim length, tc: control delay
% sched options:
% 'RRmeas' -- round-robin for meas, with 1 unused slot (match Ts of MX)
% 'MX_noACK', 'IL_noACK' 
% 'MX_piggyback', 'IL_piggyback' -- ACK sched matches measurements
% 'MX_ACK1', 'IL_ACK1' -- ACK 1 slot after control (assumes tc=1)
%
% 'SISO_' options for _piggyback or _noACK:
% 'SISOALL' - [1],[1] (std. discrete-time)
% 'SISO2' - [1 0], [0 1] for Pi_c, Pi_m
% 'SISO4' - [1 0 0 0], [0 0 1 0] for Pi_c, Pi_m
%
% 'tac' = Nv x 1 vector of time between planned control RX and ACK TX

% BR, 4/23/2014

% lots of hardcoding...  make more modular?

if(strfind(sched,'SISOALL'))
    Ts = 1;
    Pi_cBase = 1;
    Pi_mBase = 1;
    if(strfind(sched,'piggyback'))
        Pi_aBase = Pi_mBase;
    elseif(strfind(sched,'noACK'))
        Pi_aBase = 0;
    end
end

if(strfind(sched,'SISO2'))
    Ts = 2;
    Pi_cBase = [1 0];
    Pi_mBase = [0 1];
    if(strfind(sched,'piggyback'))
        Pi_aBase = Pi_mBase;
    elseif(strfind(sched,'noACK'))
        Pi_aBase = [0 0];
    end
end

if(strfind(sched,'SISO2ALLCONTROL'))
    Ts = 2;
    Pi_cBase = [1 1];
    Pi_mBase = [0 1];
    if(strfind(sched,'piggyback'))
        Pi_aBase = Pi_mBase;
    elseif(strfind(sched,'noACK'))
        Pi_aBase = [0 0];
    end
end

if(strfind(sched,'SISO4'))
    Ts = 4;
    Pi_cBase = [1 0 0 0];
    Pi_mBase = [0 0 1 0];
    if(strfind(sched,'piggyback'))
        Pi_aBase = Pi_mBase;
    elseif(strfind(sched,'noACK'))
        Pi_aBase = [0 0 0 0];
    end
end


if(strfind(sched,'RRmeas'))
    sched = 'RRmeas_noACK';
    Ts = Nv+1;  % set RR sched equal to MX sched length
    Pi_cBase = ones(Nv,Ts);
    Pi_mBase = zeros(Nv,Ts);
    for i = 1:Nv
        Pi_mBase(i,i) = 1;
    end
    Pi_aBase = zeros(Nv,Ts);
end

if( ~isempty((strfind(sched,'noACK'))) ||...
        ~isempty((strfind(sched,'piggyback'))) )

    if(strfind(sched,'block'))

        Ts = Nv + Nv;
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,Nv+i) = 1;
        end
        Pi_mBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_mBase(i,i) = 1;
        end
        
    elseif(strfind(sched,'IL'))

        Ts = Nv + Nv;
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,(2*i)) = 1;
        end
        Pi_mBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_mBase(i,(2*i-1)) = 1;
        end
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + 1;
        Pi_cBase = zeros(Nv,Ts);
        Pi_cBase(:,Ts) = ones(Nv,1);
        Pi_mBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_mBase(i,i) = 1;
        end
        
    end
    
    % set Pi_a
    if(strfind(sched,'piggyback'))
        Pi_aBase = Pi_mBase;
    else
        Pi_aBase = zeros(size(Pi_mBase));
    end

elseif(strfind(sched,'ACK1'))   % dedicated ACK slot right after u
    
    if(strfind(sched,'IL'))
        Ts = Nv + Nv + Nv;
        Pi_aBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_aBase(i,(3*i)) = 1;
        end        
        Pi_cBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_cBase(i,(3*i-1)) = 1;
        end
        Pi_mBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_mBase(i,(3*i-2)) = 1;
        end        
        
    elseif(strfind(sched,'MX'))
        Ts = Nv + Nv + 1;
        
        Pi_cBase = zeros(Nv,Ts);
        Pi_cBase(:,Nv+1) = ones(Nv,1);
        Pi_mBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_mBase(i,i) = 1;
        end        
        Pi_aBase = zeros(Nv,Ts);
        for i = 1:Nv
            Pi_aBase(i,Nv+1+i) = 1;
        end
        
    end
    
end



% compute tap
tac = zeros(Nv,1);
if(isempty(strfind(sched,'noACK')))
    for i = 1:Nv
        Pi_cPos = find(Pi_cBase(i,:))+tc;   % planned control RX
        Pi_aPos = find(Pi_aBase(i,:));
        if(Pi_cPos>Pi_aPos)
            Pi_aPos = Pi_aPos + Ts;
        end
        try
            tac(i) = Pi_aPos - Pi_cPos;
        catch
            % (SISO2ALLCONTROL -- fix eventually)
            tac(i) = 0;
        end
    end
end

nPers = ceil(Ns/Ts);
Pi_c = repmat(Pi_cBase,1,nPers);
Pi_m = repmat(Pi_mBase,1,nPers);
Pi_a = repmat(Pi_aBase,1,nPers);
Pi_c = Pi_c(:,1:Ns);
Pi_m = Pi_m(:,1:Ns);
Pi_a = Pi_a(:,1:Ns);

