function L = postOGPW1(theta,X,Y,Xd,c,c1,C1,c2,C2,C12)

Yh = heart_f(X,theta);
[Yhd, dYhd,~,dYhdx] = heart_f(Xd,theta);

if sum(imag(Yhd)) > 10^(-9)
    L = -inf;
else
    h =  c1*dYhd'+c2*dYhdx';
    H = (dYhd*C1*dYhd' + dYhd*C12*dYhdx'+ (dYhd*C12*dYhdx')'+(dYhdx*C2*dYhdx')+10^(-6)*eye(3));
    
    ROGP  = c -  h*H^(-1)*(h)';
    ROGP = ROGP/mean(diag(ROGP));
    
    SOGP = 0.02^2*ROGP+0.001^2*eye(size(X,1));
    [C,p] = chol(.5*(SOGP+SOGP'));
    
    Ldet = 2*sum(log(diag(C)));
    num = (Y-Yh)'*SOGP^(-1)*(Y-Yh);
    
    L = -1/2*Ldet-1/2*num;
end
end

