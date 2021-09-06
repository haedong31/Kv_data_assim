function  S = corr_matern3(phi,X1,X2)
%     Multidemenional separable matern correlation function with nu = 3/2
%     and lengthscale parameter phi
%
%     S = corr_matern3(phi,X1,X2) returns an covaraince matrix between the
%     sets of d demensional points X1 and X2. If length(phi) = d, then the
%     model is anisotropic. Else, the the model is anisotropic with
%     lengthcal parameter phi(1)
%
%
% phi   :  lengthscale parameters in the correlation function X1    :  n1*d
% matrix of locations X2    :  n2*d matrix of locations S     :  Covariance
% matrix



n1 = size(X1,1);
d = size(X1,2);
n2 = size(X2,1);

if length(phi) ~= d
    phi = phi(1)*ones(d,1);
end


X1_scaled = X1./repmat(phi,1,n1)';
X2_scaled = X2./repmat(phi,1,n2)';
F(1,:,:) = X2_scaled';
diff_val = abs(repmat(X1_scaled,[1,1,n2])-repmat(F,[n1,1,1]));

S = squeeze(prod(1+diff_val,2)).*exp(-squeeze(sum(diff_val,2)));

S = reshape(S,n1,n2);
