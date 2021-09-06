function sample = MH_me(U,G,L, x,M)
post_val = U(x);
sample = ones(M,length(x));
for k = 1:M
for i = 1:L
    x_p =x+G*normrnd(0,1,size(x));
    post_val_p = U(x_p);
    a = min(1,exp(-post_val_p+post_val));
    r = rand();
    if a>r || isnan(post_val_p)
        x = x_p;
        post_val = post_val_p;
        %disp('accept')        
    else
       % disp('reject')
    end
end
sample(k,:)=x;
end
end