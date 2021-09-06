function  S = corr_ddmatern3(phi,X1,X2,varargin)
%     Multidemenional separable matern correlation function with nu = 3../2
%     and lengthscale parameter phi
%
%     S = corr_matern3(phi,X1,X2) returns an covaraince matriX1 between the
%     sets of d demensional points X1 and X2. If length(phi) = d, then the
%     model is anisotropic. Else, the the model is anisotropic with
%     lengthcal parameter phi(1)
%
%
% phi   :  lengthscale parameters in the correlation function X1    :  n1....*d
% matriX1 of locations X2    :  n2...*d matriX1 of locations S     :  Covariance
% matriX1

if length(varargin) == 1
    k = varargin{1};
    Xl1 = X1(:,k);
    Xl2 = X2(:,k);
else
    k = 1;
end
n1 = size(X1,1);
n2 = size(X2,1);

if size(X1,2) > 1.5
    n1 = size(X1,1);
    n2 = size(X2,1);
    
    indo = [1:(k-1) (k+1):size(X1,2)];
    X1 = X1(:,indo);
    X2 = X2(:,indo);
    
    X1_scaled = X1./repmat(phi(indo),1,n1)';
    X2_scaled = X2./repmat(phi(indo),1,n2)';
    F(1,:,:) = X2_scaled';
    diff_val = abs(repmat(X1_scaled,[1,1,n2])-repmat(F,[n1,1,1]));
    
    S = squeeze(prod(1+diff_val,2)).*exp(-squeeze(sum(diff_val,2)));
    S = reshape(S,n1,n2);
    
    Xl1 = repmat(Xl1,[1,n2]);
    Xl2 = repmat(Xl2,[1,n1])';    
else
    S = ones(n1,n2);
    Xl1 = X1(:,k);
    Xl2 = X2(:,k);
end
Xl1 = repmat(Xl1,[1,n2]);
Xl2 = repmat(Xl2,[1,n1])';

S = S.*(exp(-abs(Xl1 - Xl2)/phi(k)).*(phi(k) -1*abs(Xl2-Xl1)))/phi(k)^3;
