function [Xh,P] = JLSKF(Xh,P,y,Uin,Dh,Nx,Nv,Nu,Np,S,A,Bu,E1,M,C,W,V,...
    ut,alphaBar,covPriorAdj)
% (uTilde1,uTilde2) (or have flexible sized uTilde matrix?)
% (tNoACK)


% covariance prior (standard)
Ppre = A*P*A' + W ; %P_{t+1|t}

if(covPriorAdj)

% scalar systems:  multiple steps
% (need to load-in "lookup table")
% work on naming vs. paper

% NOTE -- uTilde1 and uTilde2 are from prev. COMPUTED buffers
% they are not necessarily shifts of the buffer estimate in Xh
dU1 = (ut - uTilde1)*(ut - uTilde1)'; 
dU2 = (ut - uTilde2)*(ut - uTilde2)'; 

if(tNoACK==1)
    Pstar1 = (-ab*(ab-1))*dU1;
    Pstar2 = 0;
elseif(tNoACK>1)
    if(tNoACK>2)
        disp('ERROR - ACKDropped too large, using Pstar2')
    end
    Pstar1 = (- ab^4 + 2*ab^3 - 2*ab^2 + ab)*dU1^2;
    Pstar2 = (- ab^4 + 4*ab^3 - 5*ab^2 + 2*ab)*dU2^2;
else
    % no additional cov. 
    Pstar1 = 0; Pstar2 = 0;
end
Ppre = Ppre0+Pstar1+Pstar2;




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
