function L = postOGPL2(theta,X,Y,Xd,c,c1,C1)

Yh = heart_f(X,theta);
[Yhd, dYhd] = heart_f(Xd,theta);

if sum(imag(Yhd)) > 10^(-9)
    L = -inf;
else
    ROGP  = c -  c1*dYhd'*(dYhd*C1*dYhd'+10^(-6)*eye(3))^(-1)*(c1*dYhd')';%
    ROGP = ROGP/mean(diag(ROGP));
    
    SOGP = 0.02^2*ROGP+0.001^2*eye(size(X,1));
    [C,p] = chol(.5*(SOGP+SOGP'));
    
    Ldet = 2*sum(log(diag(C)));
    num = (Y-Yh)'*SOGP^(-1)*(Y-Yh);
    
    L = -1/2*Ldet-1/2*num;
end
end

