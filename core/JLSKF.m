function [Xh,P] = JLSKF(Xh,P,y,Uin,Dh,Nx,Nv,Nu,Np,S,A,Bu,E1,M,C,W,V,...
    uOptions,ab,covPriorAdj,tNoACK,t,tp)

% add more Pstars...
% more streamlined approach vs. hardcode each # steps?

% covariance prior (standard)
Ppre0 = A*P*A' + W ; %P_{t+1|t}

if( covPriorAdj && (tNoACK>0) )
    
    if(size(A,1)==1)
        % scalar systems:  multiple steps
        % (need to load-in "lookup table")
        % work on naming vs. paper
        if(tNoACK>2)
            fprintf('\nt=%d, tp=%d, KF ERROR - ACKDropped too large, using Pstar2\n',t,tp)
            tNoACK=2;
        end
        
        % NOTE -- uHistory are from prev. COMPUTED buffers
        % they are not necessarily shifts of the buffer estimate in Xh
        
        % SENT control command for step t
        ut = E1*Uin;
        
        dU = zeros(1,tNoACK);
        for i = 1:tNoACK
            dU(:,i) = (ut - uOptions(:,i))*(ut - uOptions(:,i))';
        end
        
        Pstar = zeros(1,10);
        if(tNoACK==1)
            
            Pstar(1) = (-ab*(ab-1))*dU(:,1);
            %Pstar2 = 0;
            fprintf('\nt=%d, tp=%d, Pstar1 = %f \n',t, tp, Pstar(1))
            
        elseif(tNoACK==2)
            
            % uses diff with 1-step prev. plan
            Pstar(1) = (- ab^4 + 2*ab^3 - 2*ab^2 + ab)*dU(:,2)^2;
            % uses diff with 2-step prev. plan (earliest)
            Pstar(2) = (- ab^4 + 4*ab^3 - 5*ab^2 + 2*ab)*dU(:,1)^2;
            fprintf('\nt=%d, tp=%d, Pstar1 = %f, Pstar2 = %f \n',t,tp,Pstar(1),Pstar(2))

        elseif(tNoACK>2)
            
            % (fill in with more Pstars...)
            %Pstar1 = 0.00111; Pstar2 = 0.00222; Pstar3 = 0.00333; 
            disp('(Pstar>2)')
            
        end
        
        Ppre = Ppre0+sum(Pstar);
        
    else
        disp('WARNING - cov. prior adjust not formulated for multivar. systems')
        
        %{
        % single-step: Pstar
        % works for MIMO
        % covariance prior: scale control by alphabar -- diagonal matrix
        ubuff = E1*M*Xh(Nx+1:end);    % one shift of previous buffer
        qq = (ut - ubuff)*(ut - ubuff)' ;
        ab = diag(alphaBar);
        Z = ab*qq*ab ;   % off diagonal elements
        dum = ab*qq ;
        for i = 1:Nu; Z(i,i) = dum(i,i) ; end;
        Padd =  Bu*(Z + qq - ab*qq - qq*ab)*Bu';
        Ppre = Ppre + Padd;
        %}
        
    end
else
    
    Ppre = Ppre0;
    
end

% Kalman gain = fcn of Ppre
L = ((Ppre*C')/(C*Ppre*C'+V))*S;   % L_{t}
P = Ppre - L*(C*Ppre);

% propagate system
I = eye(size(A));
AAHat = [(I-L*C)*A,(I-L*C)*Bu*E1*M*(eye(Np*Nu*Nv)-Dh);...
    zeros(Nv*Np*Nu,Nx),M*(eye(Np*Nu*Nv)-Dh)];
BUHat = [(I-L*C)*Bu*E1*Dh;Dh];
Xh = AAHat*Xh + BUHat*Uin+[L;zeros(Np*Nu*Nv,size(y,1))]*y;

end
