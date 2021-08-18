function [yksum, comp_currents] = kcurrent_model(model_struct)
    current_names = fieldnames(model_struct);
    num_currents = length(current_names);
    
    tune_params = cell(num_currents, 1);
    for i = 1:num_currents
        tune_params = model_struct.(current_names{i});
    end


end

function [current_trace] = gen_matching_current(current_name, tune_param)
    switch current_name
        case "ikto"
            % generate ikto
            num_kto_param = 17;
            fixed_kto_idx = setdiff(1:num_kto_param, tune_param);

            
        case "ikslow1"
            % generate ikslow1
        case "ikslow2"
            % generate ikslow2

        case "ikur"
            % generate ikur
        case "ikss"
            % generate ikss
            ff
        case "ik1"
            % generate ik1
            ff
    end
end
