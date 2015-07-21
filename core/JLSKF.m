function [XHat_out,P_out,Pstar_out] = JLSKF(XHat_in,P_in,y_KF,U_sent,...
    Dc_hat_all,N_STATES,N_AGENTS,N_CONTROLS,N_HORIZON,Dm_now,...
    A,Bu,E,M,C,W,V,ALPHAC_BAR,cov,pd)
% Kalman filter for missed measurements and controls
% [XHat_out,P_out] = JLSKF(XHat_in,P_in,y_KF,U_sent,Dc_hat_all,...
%     N_STATES,N_AGENTS,N_CONTROLS,N_HORIZON,Dm_now,A,Bu,E,M,C,W,V,...
%     ALPHAC_BAR,cov,pd)
% XHat is state estimate (system + buffer), all time steps
% P is error cov., all time steps
% y_KF: meas into KF
% Dm_now: meas success matrix
% U_sent: control plan as sent (used if packet successful)
% Dc_hat: estimate of control jump variable
% W, V: process, meas. noise covariances
% ALPHAC_BAR: control packet success probability
% cov=v struct:  (if cov==0 or [], cov_prior_adjust is not used)
% else, cov must be a struct with the following fields:
%   u_options: N_CONTROLS x t_NoACK vector of alternate control actions 
%       u_options is indexed backwards
%       u_options(:,1) is most recent, (:,t_NoACK_KF) is furthest back
%   t_NoACK: time since ACK received (for cov_prior_adj)
%   Pstar_coefficients: vector of coefficients for Pstar(1:n_stars-1)
%   Pstar_final_coefficients: vector of last term coeffs: Pstar(n_stars)
% pd struct: if(pd==0 or [], print_debug_KF = 0)
%   t: true time at estimator when fcn is called (mostly for debugging)
%   t_KF: physical time step KF is updating -- x_hat(t_KF|t_KF) posterior

if(isstruct(cov))
    
    cov_prior_adj = 1;
    u_options = cov.u_options;
    t_NoACK = cov.t_NoACK;
    if(size(Bu,2)==1)
        Pstar_coefficients = cov.Pstar_coefficients;
        Pstar_final_coefficients = cov.Pstar_final_coefficients;
        n_stars = length(Pstar_coefficients);
    end
elseif(isempty(cov) || cov==0)
    
    cov_prior_adj = 0;
    
end

Pstar_out = [];

if( isstruct(pd) )
    
    t = pd.t; t_KF = pd.t_KF;
    print_debug_kf = 1;
    
elseif(isempty(pd) || pd==0)
    
    print_debug_kf = 0;
    
end

% P_{t+1|t} covariance prior (standard)
% note - if needed, this W should be scaled by Bw (as in setupSystemJLSPPC) 
P_pre_0 = A*P_in*A' + W ; 

if( cov_prior_adj && (max(t_NoACK)>0) )
    
    % SENT control command for step t
    u_t = E*U_sent;
    
    q = zeros(size(u_t,1),max(t_NoACK));
    for i = 1:max(t_NoACK)
        q(:,i) = (u_t - u_options(:,i));
    end
    
    % Single input systems
    if(size(Bu,2)==1)
        
        % tNoACK is a scalar
        if(length(t_NoACK)>1)
            disp('ERROR: tNoACK is vector, expected scalar')
        end
        
        if(t_NoACK>n_stars)
            if(print_debug_kf)
                fprintf('\nt=%d, KF tKF=%d ERROR: tNoACK too large, using Pstar%d\n',...
                    t,t_KF,n_stars)
            end
            t_NoACK=n_stars;
        end
        
        % NOTE -- uHistory are from prev. COMPUTED buffers
        % they are not necessarily shifts of the buffer estimate in Xh
        Pstar = zeros(size(P_in,1),size(P_in,2),t_NoACK);
        
        % Construct covariance addition terms: P*
        if(t_NoACK>1)
            % coefficients for all except last term:
            for j = 1:(t_NoACK-1)
                Pstar(:,:,j) = Pstar_coefficients(j)*...
                    Bu*q(:,j)*q(:,j)'*Bu';
            end
        end
        % separate coefficient form for the last term:
        Pstar(:,:,t_NoACK) = Pstar_final_coefficients(t_NoACK)*...
            Bu*q(:,t_NoACK)*q(:,t_NoACK)'*Bu';
        
        if(print_debug_kf)
            
            % print out sum of Pstar terms
            
            % notes on timing debug output: 
            % t is the global time (at estimator)
            % t_KF is the a posteriori time step that the KF is updating
            % tNoACK_KF is the counter for the priors at this step 
            %   eg --   if u(t_KF-1) has been ACK'd, tNoACK_KF=0
            %           if not, tNoACK_KF>0    
            
            if(size(Pstar,1)==1)
                fprintf('\nt=%d, KF tKF=%d, tNoACK_KF = %d, P* = %f \n',...
                    t, t_KF,t_NoACK,squeeze(sum(Pstar,3)))
            else
                fprintf('\nt=%d, KF tKF=%d, tNoACK_KF = %d, P* = \n',...
                    t, t_KF,t_NoACK)
                disp(squeeze(sum(Pstar,3)))
            end
            
            % (deeper debug option of printing out individual terms)
            % disp(Pstar)

        end
        
        % add on the new terms
        P_pre = P_pre_0 + sum(Pstar,3);
        
    else
        
        % single-step: Pstar
        
        % (zero out dU manually for ACK'd channels - need to test)
        for i = 1:length(t_NoACK)
            if(t_NoACK(i) == 0)
                q(i,1) = 0;
            end
        end
        
        % off diagonal elements
        E_AZA = diag(ALPHAC_BAR)*q(:,1)*q(:,1)'*diag(ALPHAC_BAR);
        
        % replace diagonal elements
        dum = diag(ALPHAC_BAR)*q(:,1)*q(:,1)';
        for i = 1:N_CONTROLS; E_AZA(i,i) = dum(i,i) ; end;
        
        Pstar = Bu*(E_AZA - diag(ALPHAC_BAR)*q(:,1)*q(:,1)'*...
            diag(ALPHAC_BAR))*Bu';
        
        P_pre = P_pre_0 + Pstar;
        
        if(print_debug_kf)
            fprintf('\nt=%d, KF tKF=%d, MIMO P*: \n',t,t_KF)
            disp(Pstar)
        end
        
        
    end
    
    Pstar_out = Pstar;
    
else
    
    P_pre = P_pre_0;
    
end

% Kalman gain = fcn of Ppre
L = ((P_pre*C')/(C*P_pre*C'+V))*Dm_now;   % L_{t}
P_out = P_pre - L*(C*P_pre);

% propagate system
I = eye(size(A));
AAHat = [(I-L*C)*A,...
    (I-L*C)*Bu*E*M*sparse(eye(N_HORIZON*N_CONTROLS*N_AGENTS)-Dc_hat_all);...
    zeros(N_AGENTS*N_HORIZON*N_CONTROLS,N_STATES),...
    M*sparse(eye(N_HORIZON*N_CONTROLS*N_AGENTS)-Dc_hat_all)];
BUHat = [(I-L*C)*Bu*E*Dc_hat_all;Dc_hat_all];
XHat_out = AAHat*XHat_in + BUHat*U_sent + ...
    [L;zeros(N_HORIZON*N_CONTROLS*N_AGENTS,size(y_KF,1))]*y_KF;

end
