function [z] = obj_rmse_over_volt(p, hold_volt, hold_idx, volt, t, exp_ksum, Ek, param_select)
    rmse = zeros(length(volt),1);
    for i=1:length(volt)
        volt_idx=i;
        time_space = cell(1,3);
        time_space{1} = t;
        time_space{2} = t(1:hold_idx(volt_idx));
        time_space{3} = t(hold_idx(volt_idx)+1:end) - t(hold_idx(volt_idx)+1);
        rmse(i) = obj_rmse(p, hold_volt, volt(i), time_space, Ek, exp_ksum(:,i), param_select); 
    end
    z = sum(rmse);
end