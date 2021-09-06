function [ov,dov,ovdx,dovdx]  = heart_f(T,nu)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


t_1 = exp(nu(1));
t_2 = exp(nu(2));
t_3 = exp(nu(3));

A = [-t_3-t_2 t_1                       0           0;%O
    t_2  -t_1-t_2                  t_1         0; %C1     
     0           t_2                -t_1-t_2    t_1; %C2
    0               0                 t_2         -t_1]; %C3

y0 = [0 0  0 1]';

dA{1} = [1 0 0 0;
    -1 0 1 0;    
    0 0 -1 1;
    0 0 0 -1];
dA{2}= [0 -1 0 0;
    -1 1 0 0;    
    1 0 -1 0;
    0 0 1 0];
dA{3} = [0 -1 0 0;
    0 0 0 0;    
    0 0 0 0;
    0 0 0 0];
obs = [1 0 0 0];
[ov,dov,ovdx,dovdx] = dexpme2(A,dA,y0,obs,exp(T));
ovdx = exp(T).*ovdx;
dovdx = repmat(exp(T),1,length(nu)).*dovdx;

if nargout>1
    dov = (dov.*repmat(exp(nu)',length(T),1))';
    dovdx = (dovdx.*repmat(exp(nu)',length(T),1))';
end

end
%*********************************************************************************************
function [F,dF,Fdx,dFdx,dy0]=dexpme2(M,dM,y0,obs,t) % by Lubomír Bran?ík, 2008
% Eigenvalues/eigenvectors decomposition - usage of Matlab eig
[U,Q]=eig(M);
eQrS=exp(diag(Q));
Rem2 = (U\y0);
Rem3 = (U^(-1));
Rem1 = (obs*U);

Rem4 = (obs*M*U);

F = ones(length(t),1);
Fdx = ones(length(t),1);
dF = ones(length(t),length(dM));
dFdx = ones(length(t),length(dM));
dy0 = ones(length(t),length(M));

for k = 1:length(t)
    eQr = eQrS.^t(k);
    F(k) = Rem1*diag(eQr)*Rem2;
    Fdx(k) = Rem4*diag(eQr)*Rem2;
end
if nargout >= 2
    for l = 1:length(dM)
        Rem5 = (obs*dM{l}*U);
        H = U\dM{l}*U;
        for k = 1:length(t)
            eQr = eQrS.^t(k);
            P=diag(t(k)*diag(H).*eQr);
            for i=1:length(H)-1
                for j=i+1:length(H)
                    RQ=(eQr(i)-eQr(j))/(Q(i,i)-Q(j,j));
                    P(i,j)=H(i,j)*RQ; P(j,i)=H(j,i)*RQ;
                end
            end
            dF(k,l) =Rem1*P*Rem2;
            dFdx(k,l) =Rem4*P*Rem2+Rem5*diag(eQr)*Rem2;
        end
    end
end
if nargout >=3    
    for k = 1:length(t)
        eQr = eQrS.^t(k);
        dy0(k,:) = Rem1*diag(eQr)*Rem3 ;
    end
end
end
%*********************************************************************************************
