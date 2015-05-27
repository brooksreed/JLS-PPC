function [Xh,P] = JLSKF(Xh,P,y,Uin,Dh,Nx,Nv,Nu,Np,S,A,Bu,E1,M,C,W,V,...
    ut,alphaBar,covPriorAdj)

% covariance prior (standard)
Ppre = A*P*A' + W ; %P_{t+1|t}

% covariance prior: scale control by alphabar -- diagonal matrix
if(covPriorAdj)
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
