function z = obj_biomarkers(p, kcurrent_model, model_struct, volt_space, time_space, yksum)
    volts = volt_space{2};
    num_volts = length(volts);

    protocol = cell(4, 1);
    protocol{1} = volt_space{1};
    protocol{3} = time_space;
    protocol{4} = volt_space{3};
    
    mse_list = NaN(num_volts, 1);
    for i = 1:num_volts
        protocol{2} = volts(i);
        
        % experimental data
        yksum_i = yksum(:, i);
        mexp = get_biomarkers(yksum_i, time_space);

        % computer model
        [yksum_hat, ~] = kcurrent_model(p, model_struct, protocol);
        mhat = get_biomarkers(yksum_hat, time_space);

        mse_list(i) = mean((mexp - mhat).^2);
    end
    z = sum(mse_list);
end

function m = get_biomarkers(y, time_space)
    hold_idx = time_space{4};
    ysub = y((hold_idx+1):end);
    m = NaN(3, 1);

    % peak 
    yearly = ysub(1:floor(length(ysub)*0.05));
    m(1) = max(yearly);

    % middle
    m(2) = ysub(floor(length(ysub)/2));
    
    % tail
    m(3) = y(end);
end
