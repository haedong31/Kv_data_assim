function z = obj_biomarkers(p, kcurrent_model, model_struct, volt_space, time_space, yksum)
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};
    
    diff_list = NaN(num_volts, 1);
    for i = 1:num_volts
        protocol{2} = volts(2);
        
        % experimental data
        yksum_i = yksum(:, i);
        [m1_exp, m2_exp] = get_biomarkers(yksum_i, time_space);

        % computer model
        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        [m1_hat, m2_hat] = get_biomarkers(yksum_hat, time_space);

        if (isequal(m2_exp, -1) || isequal(m2_hat, -1))
            diff_list(i) = abs(m1_exp-m1_hat);
        else
            diff_list(i) = abs(m1_exp-m1_hat) + abs(m2_exp-m2_hat);
        end
    end
    z = sum(diff_list);
end

function [m1, m2] = get_biomarkers(y, time_space)
    t = time_space{1};
    hold_idx = time_space{4};

    % peak
    ysub1 = y((hold_idx+1):end);
    tsub1 = t((hold_idx+1):end);

    yearly = ysub1(1:floor(length(ysub1)*0.05));
    [m1, m1_idx] = max(yearly);
    
    % time constant
    ysub2 = ysub1(m1_idx:end);
    tsub2 = tsub1(m1_idx:end);
    
    tau_current = m1*exp(-1);
    if tau_current < y(end)
        m2 = -1;
    else
        [~, m2_idx] = min(abs(tau_current - ysub2));
        m2 = tsub2(m2_idx);
    end
end
