function L = postKO(theta,X,Y,RKO)
Yh = heart_f(X,theta);

if sum(imag(Yh)) > 10^(-9)
    L = -inf;
else
    SKO = 0.02^2*RKO+0.001^2*eye(size(X,1));
    [C,p] = chol(.5*(SKO+SKO'));
    
    Ldet = 2*sum(log(diag(C)));
    num = (Y-Yh)'*SKO^(-1)*(Y-Yh);
    
    L = -1/2*Ldet-1/2*num;
end
end

