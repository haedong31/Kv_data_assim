function z = obj_rmse(p, kcurrent_model, mdl_struct, pdefault, protocol, yksum)
    volts = protocol{3};
    num_volts = length(volts);
    hold_idx = length(protocol{5});

    rmse_list = NaN(num_volts,1);
    for i = 1:num_volts
        yi = yksum(:,i);
        ymx = kcurrent_model(p, mdl_struct, pdefault, protocol, volts(i));
        yhat = sum(ymx,2);
        rmse_list(i) = sqrt(mean((yi(hold_idx+1:end) - yhat(hold_idx+1:end)).^2));
    end
    z = sum(rmse_list);
end
