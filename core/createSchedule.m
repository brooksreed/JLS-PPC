function [PI_C,PI_M,PI_A,TAU_AC,T_S] = createSchedule(sched,N_VEH,...
    SIM_LEN,TAU_C)
% creates schedule time series variables and TAU_AC and T_S parameters
% PI_C,PI_M,PI_A,TAU_AC,T_S] = createSchedule(sched,N_VEH,SIM_LEN,TAU_C)
%
% TAU_AC = N_VEH x 1 vector of time between planned control RX and ACK TX
% T_S: schedule length (one period)
%
% sched options:
% 'RRmeas' -- round-robin for meas, with 1 unused slot (match Ts of MX)
% 'MX_noACK', 'IL_noACK' 
% 'MX_piggyback', 'IL_piggyback' -- ACK sched matches measurements
% 'MX_ACK1', 'IL_ACK1' -- ACK 1 slot after control (assumes tc=1)
%
% 'SISO_' options for _piggyback or _noACK:
% 'SISOALL' - [1],[1] (std. discrete-time)
% 'SISO2' - [1 0], [0 1] for PI_C, PI_M
% 'SISO4' - [1 0 0 0], [0 0 1 0] for PI_C, PI_M
%

% TO DO: 
% lots of hardcoding...  make more modular?

if(strfind(sched,'SISOALL'))
    T_S = 1;
    Pi_c_base = 1;
    Pi_m_base = 1;
    if(strfind(sched,'piggyback'))
        Pi_a_base = Pi_m_base;
    elseif(strfind(sched,'noACK'))
        Pi_a_base = 0;
    end
end

if(strfind(sched,'SISO2'))
    T_S = 2;
    Pi_c_base = [1 0];
    Pi_m_base = [0 1];
    if(strfind(sched,'piggyback'))
        Pi_a_base = Pi_m_base;
    elseif(strfind(sched,'noACK'))
        Pi_a_base = [0 0];
    end
end

if(strfind(sched,'SISO2ALLCONTROL'))
    T_S = 2;
    Pi_c_base = [1 1];
    Pi_m_base = [0 1];
    if(strfind(sched,'piggyback'))
        Pi_a_base = Pi_m_base;
    elseif(strfind(sched,'noACK'))
        Pi_a_base = [0 0];
    end
end

if(strfind(sched,'SISO4'))
    T_S = 4;
    Pi_c_base = [1 0 0 0];
    Pi_m_base = [0 0 1 0];
    if(strfind(sched,'piggyback'))
        Pi_a_base = Pi_m_base;
    elseif(strfind(sched,'noACK'))
        Pi_a_base = [0 0 0 0];
    end
end


if(strfind(sched,'RRmeas'))
    sched = 'RRmeas_noACK';
    T_S = N_VEH+1;  % set RR sched equal to MX sched length
    Pi_c_base = ones(N_VEH,T_S);
    Pi_m_base = zeros(N_VEH,T_S);
    for i = 1:N_VEH
        Pi_m_base(i,i) = 1;
    end
    Pi_a_base = zeros(N_VEH,T_S);
end

if( ~isempty((strfind(sched,'noACK'))) ||...
        ~isempty((strfind(sched,'piggyback'))) )

    if(strfind(sched,'block'))

        T_S = N_VEH + N_VEH;
        Pi_c_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_c_base(i,N_VEH+i) = 1;
        end
        Pi_m_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_m_base(i,i) = 1;
        end
        
    elseif(strfind(sched,'IL'))

        T_S = N_VEH + N_VEH;
        Pi_c_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_c_base(i,(2*i)) = 1;
        end
        Pi_m_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_m_base(i,(2*i-1)) = 1;
        end
        
    elseif(strfind(sched,'MX'))
        T_S = N_VEH + 1;
        Pi_c_base = zeros(N_VEH,T_S);
        Pi_c_base(:,T_S) = ones(N_VEH,1);
        Pi_m_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_m_base(i,i) = 1;
        end
        
    end
    
    % set Pi_a
    if(strfind(sched,'piggyback'))
        Pi_a_base = Pi_m_base;
    else
        Pi_a_base = zeros(size(Pi_m_base));
    end

elseif(strfind(sched,'ACK1'))   % dedicated ACK slot right after u
    
    if(strfind(sched,'IL'))
        T_S = N_VEH + N_VEH + N_VEH;
        Pi_a_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_a_base(i,(3*i)) = 1;
        end        
        Pi_c_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_c_base(i,(3*i-1)) = 1;
        end
        Pi_m_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_m_base(i,(3*i-2)) = 1;
        end        
        
    elseif(strfind(sched,'MX'))
        T_S = N_VEH + N_VEH + 1;
        
        Pi_c_base = zeros(N_VEH,T_S);
        Pi_c_base(:,N_VEH+1) = ones(N_VEH,1);
        Pi_m_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_m_base(i,i) = 1;
        end        
        Pi_a_base = zeros(N_VEH,T_S);
        for i = 1:N_VEH
            Pi_a_base(i,N_VEH+1+i) = 1;
        end
        
    end
    
end



% compute tap
TAU_AC = zeros(N_VEH,1);
if(isempty(strfind(sched,'noACK')))
    for i = 1:N_VEH
        Pi_c_positive = find(Pi_c_base(i,:)) + TAU_C;   
        Pi_a_ = find(Pi_a_base(i,:));
        if(Pi_c_positive>Pi_aPos)
            Pi_aPos = Pi_aPos + T_S;
        end
        try
            TAU_AC(i) = Pi_aPos - Pi_c_positive;
        catch
            % (SISO2ALLCONTROL -- fix eventually)
            TAU_AC(i) = 0;
        end
    end
end

n_periods = ceil(SIM_LEN/T_S);
PI_C = repmat(Pi_c_base,1,n_periods);
PI_M = repmat(Pi_m_base,1,n_periods);
PI_A = repmat(Pi_a_base,1,n_periods);
PI_C = PI_C(:,1:SIM_LEN);
PI_M = PI_M(:,1:SIM_LEN);
PI_A = PI_A(:,1:SIM_LEN);

