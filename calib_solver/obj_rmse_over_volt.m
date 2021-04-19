function [z] = obj_rmse_over_volt(p, hold_volt, volt, time_space, Ek, exp_ksum)
    rmse = zeros(length(volt),1);
    for i=1:length(volt)
       rmse(i) = obj_rmse(p, hold_volt, volt(i), time_space, Ek, exp_ksum(:,i)); 
    end
    z = sum(rmse);
end