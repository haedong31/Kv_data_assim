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
        mexp = get_biomarkers(yksum_i, time_space);

        % computer model
        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        mhat = get_biomarkers(yksum_hat, time_space);

        diff_list(i) = abs(mexp-mhat);
    end
    z = sum(diff_list);
end

function m = get_biomarkers(y, time_space)
    hold_idx = time_space{4};

    % peak
    ysub1 = y((hold_idx+1):end);
    yearly = ysub1(1:floor(length(ysub1)*0.05));
    m = max(yearly);    
end
