function [z] = obj_pretrn_rmse(p, phase, kcurrent_model, model_struct, volt_space, time_space, yksum, pretrn_mdl)
    hold_idx = time_space{4};
    end_idx = time_space{5};
    
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};

    break_pt = floor((end_idx - hold_idx)*0.05) + hold_idx;
    rmse_list = zeros(num_volts, 1);
    switch phase
        case 'early'
            for i = 1:num_volts
                yksum_i = yksum(:, i);
                protocol{2} = volts(i);
        
                [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol, pretrn_mdl);
                
                yksum_i = yksum_i(1:break_pt);
                yksum_hat = yksum_hat(1:break_pt);

                rmse_list(i) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
            end
        case 'tail'
            for i = 1:num_volts
                yksum_i = yksum(:, i);
                protocol{2} = volts(i);
        
                [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol, pretrn_mdl);

                yksum_i = yksum_i((break_pt + 1):end);
                yksum_hat = yksum_hat((break_pt + 1):end);

                rmse_list(i) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
            end
        otherwise
            for i = 1:num_volts
                yksum_i = yksum(:, i);
                protocol{2} = volts(i);
        
                [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol, pretrn_mdl);
                rmse_list(i) = sqrt(mean((yksum_i((hold_idx + 1):end) - yksum_hat((hold_idx + 1):end)).^2));
            end
    end
    z = sum(rmse_list);
end
