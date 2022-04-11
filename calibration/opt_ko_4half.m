clc
close all
clearvars
warning('off', 'all')

%----- Code arguments & Model information %-----
group = "ko";
exp_num = "exp41";
save_dir = strcat("calib_",exp_num,"_",group);

% currents to be included
current_names = {'iktof', 'ikslow1', 'ikslow2', 'ikss'};

tune_idx = cell(5,1); % iktof, iktos, ikslow1, ikslow2, ikss
tune_idx{1} = [1, 2, 4, 5, 7, 11, 13];
% tune_idx{2} = [2, 3];
tune_idx{3} = [1, 2, 4, 5, 9, 10, 11];
tune_idx{4} = [2, 3];
tune_idx{5} = [1, 2, 3, 4];

pdefault = cell(5,1);
pdefault{1} = [33, 15.5, 20, 7, 0.03577, 0.06237, 0.18064, 0.3956, ...
    0.000152, 0.067083, 0.00095, 0.051335, 0.3846];
pdefault{2} = [-1050, 270, 0.0629];
pdefault{3} = [22.5, 45.2, 40.0, 7.7, 5.7, 0.0629, 6.1, 18, 2.058, 803.0, 0.16];
pdefault{4} = [4912, 5334, 0.16];
pdefault{5} = [0.0862, 1235.5, 13.17, 0.0611];

[mdl_struct, psize] = gen_mdl_struct(current_names, tune_idx);
[p0, lb, ub] = gen_param_bounds(mdl_struct, psize, pdefault);

%----- Protocol -----%
% voltage
holdv = -70;
minv = 50;
vstep = 10;
num_steps = 1;
ek = -91.1;

% time
ideal_hold_time = 470;
ideal_end_time = 25*1000;

% save time information later when experimental data imlported
protocol = cell(6,1);
protocol{1} = holdv;
protocol{2} = ek;
protocol{3} = minv:vstep:(minv+vstep*(num_steps-1));

% meta data of experimental datasets
matching_table = readtable(fullfile(pwd,"mgat1ko_data",strcat("matching-table-",group,".xlsx")));
file_names = matching_table.trace_file_name_4half;
data_dir = fullfile("mgat1ko_data",strcat(group,"-preprocessed-25s"));

% exclude null rows
loop_idx = [];
for i = 1:size(file_names,1)
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end
num_files = length(loop_idx);

%----- Optimization loop -----%
num_iters = 5;
options = optimoptions(@fmincon, ...
    'Algorithm','interior-point', 'Display','off', ...
    'MaxFunctionEvaluations',1e+6, ...
    'SpecifyObjectiveGradient',false);
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

outf = fopen(strcat(exp_num,"_",group,".txt"), 'w');
for i = 1:num_files
    % read data
    file_path = fullfile(pwd,data_dir,file_names{loop_idx(i)});
    trace_data = table2array(readtable(file_path));
    t = trace_data(:, 1);
    yksum = trace_data(:,2:end);

    % estimate the critical time points
    [~,ideal_hold_idx] = min(abs(t - ideal_hold_time));
    [~,ideal_end_idx] = min(abs(t - ideal_end_time));
    t = t(1:ideal_end_idx);
    yksum = yksum(1:ideal_end_idx,:);

    % save time information in the protocol
    protocol{4} = t;
    protocol{5} = t(1:ideal_hold_idx);
    pulset = t(ideal_hold_idx+1:end);
    protocol{6} = pulset - pulset(1);
    
    % objective function
    opt_fun = @(p) obj_rmse(p, @kcurrent_model2, mdl_struct, pdefault, protocol, yksum);

    % run optimization
    rmse_list = zeros(num_iters, 1);
    sol_list = cell(num_iters, 1);

    % first run with p0
    try
    [sol, fval] = fmincon(opt_fun, p0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    sol_list{1} = sol;
    rmse_list(1) = fval;

    outs = sprintf("[File %i/%i] %s [Reps %i/%i] Min RMSE: %f", l, len_loop_idx, file_names{i}, 1, num_iters, fval);
    fprintf(outf, '%s\n', outs);
    disp(outs)

    % random intialization
    running_p0 = lhsdesign(num_iters, psize);
    running_p0 = scale_param(running_p0, lb, ub);
    
    for j = 2:num_iters
        % optimization
        try
            [sol, fval] = fmincon(opt_fun, running_p0(j,:), A, b, Aeq, beq, lb, ub, nonlcon, options);
            sol_list{j} = sol;
            rmse_list(j) = fval;
        catch
            rmse_list(j) = 1e+3;
        end
        outs = sprintf('[File %i/%i] %s [Reps %i/%i] Min RMSE: %f', l, len_loop_idx, file_names{i}, j, num_iters, fval);
        fprintf(outf, '%s\n', outs);
        disp(outs)
    end

    [~, best_fit_idx] = min(rmse_list);
    best_sol = sol_list{best_fit_idx};

    % save calibrated solution
    save_path = fullfile(pwd, save_dir, file_names{i});
    sol_mx = zeros(max_param_len, num_currents);
    for j = 1:num_currents
        switch current_names{j}
        case 'iktof'
            sol_kto = ktof_default;
            sol_kto(tune_idx1_ktof) = best_sol(idx_info2{j});
            sol_kto = [sol_kto, NaN(1, (max_param_len-length(sol_kto)))];
            sol_mx(:, j) = sol_kto;
        case 'ikslow1'
            sol_kslow1 = kslow1_default;
            sol_kslow1(tune_idx1_kslow1) = best_sol(idx_info2{j});
            sol_kslow1 = [sol_kslow1, NaN(1, (max_param_len-length(sol_kslow1)))];
            sol_mx(:, j) = sol_kslow1;
        case 'ikslow2'
            sol_kslow2 = kslow2_default;
            sol_kslow2(tune_idx1_kslow2) = best_sol(idx_info2{j});
            sol_kslow2 = [sol_kslow2, NaN(1, (max_param_len-length(sol_kslow2)))];
            sol_mx(:, j) = sol_kslow2;
        case 'ikss'
            sol_kss = kss_default;
            sol_kss(tune_idx1_kss) = best_sol(idx_info2{j});
            sol_kss = [sol_kss, NaN(1, (max_param_len-length(sol_kss)))];
            sol_mx(:, j) = sol_kss;
        case 'ikur'
            sol_kur = kur_default;
            sol_kur(tune_idx1_kur) = best_sol(idx_info2{j});
            sol_kur = [sol_kur, NaN(1, (max_param_len-length(sol_kur)))];
            sol_mx(:, j) = sol_kur;
        case 'ik1'
            sol_k1 = k1_default;
            sol_k1(tune_idx1_k1) = best_sol(idx_info2{j});
            sol_k1 = [sol_k1, NaN(1, (max_param_len-length(sol_k1)))];
            sol_mx(:, j) = sol_k1;
        end
    end
%     obj_rmse(best_sol, 'all', @kcurrent_model1, model_struct, volt_space, time_space, yksum)
    writematrix(string(current_names) , save_path, "Sheet","Parameters", "Range","A1");
    writematrix(sol_mx, save_path, "Sheet","Parameters", "Range","A2");
end
fclose(outf);

%% ----- Custom functions ----- %%
function [p0, lb, ub] = gen_param_bounds(mdl_struct, psize, pdefault)
    p0 = NaN(psize,1);
    lb = NaN(psize,1);
    ub = NaN(psize,1);    

    % lower bound
    lb_ktof = NaN(1,length(pdefault{1}));
    lb_ktof([1:4, 13]) = [-70, eps, eps, 1, eps];
    lb_ktof([5,6,10,12]) = pdefault{1}([5,6,10,12])*0.15;
    lb_ktof([7,8]) = pdefault{1}([7,8])*0.1;
    lb_ktof([9,11]) = pdefault{1}([9,11])*0.7;
    lb_ktos = pdefault{2};
    lb_ktos(1) = lb_ktos(1)*1.95;
    lb_ktos(2:end) = lb_ktos(2:end)*0.05;
    lb_kslow1 = [-70, -70, -70, 1, 1, eps, eps, eps, eps, 50+eps, eps];
    lb_kslow2 = [eps, 5000, eps];
    lb_kss = [eps, eps, eps, eps];

    % upper bound
    ub_ktof = NaN(1,length(pdefault{1}));
    ub_ktof([1:4, 13]) = [70, 40, 40, 14, 1];
    ub_ktof(5:12) = pdefault{1}(5:12)*1.95;
    ub_ktos = pdefault{2};
    ub_ktos(1) = ub_ktos(1)*0.05;
    ub_ktos(2:end) = ub_ktos(2:end)*1.95;
    ub_kslow1 = [50, 50, 50, 20, 20, 1, 30, 50, 20, 2000, 1];
    ub_kslow2 = [5000, 10000, 1];
    ub_kss = [1, 2000, 100, 1];

    % assign the bounds according to the model structure
    current_names = {mdl_struct.name};
    idx_info1 = {mdl_struct.idx1};
    idx_info2 = {mdl_struct.idx2};
    for i = 1:length(current_names)
        switch current_names{i}
            case 'iktof'
                p0(idx_info2{i}) = pdefault{1}(idx_info1{i});
                lb(idx_info2{i}) = lb_ktof(idx_info1{i});
                ub(idx_info2{i}) = ub_ktof(idx_info1{i});
            case 'iktos'
                p0(idx_info2{i}) = pdefault{2}(idx_info1{i});
                lb(idx_info2{i}) = lb_ktos(idx_info1{i});
                ub(idx_info2{i}) = ub_ktos(idx_info1{i});
            case 'ikslow1'
                p0(idx_info2{i}) = pdefault{3}(idx_info1{i});
                lb(idx_info2{i}) = lb_kslow1(idx_info1{i});
                ub(idx_info2{i}) = ub_kslow1(idx_info1{i});
            case 'ikslow2'
                p0(idx_info2{i}) = pdefault{4}(idx_info1{i});
                lb(idx_info2{i}) = lb_kslow2(idx_info1{i});
                ub(idx_info2{i}) = ub_kslow2(idx_info1{i});
            case 'ikss'
                p0(idx_info2{i}) = pdefault{5}(idx_info1{i});
                lb(idx_info2{i}) = lb_kss(idx_info1{i});
                ub(idx_info2{i}) = ub_kss(idx_info1{i});
        end
    end
end

function scaledp = scale_param(unitp, lb, ub)
    [num_pt, num_var] = size(unitp);
    scaledp = NaN(num_pt, num_var);
    for i = 1:num_pt
        scaledp(i, :) = unitp(i,:).*(ub-lb) + lb;
    end
end

function z = obj_rmse(p, kcurrent_model, mdl_struct, pdefault, protocol, yksum)
    volts = protocol{3};
    num_volts = length(volts);
    hold_idx = length(protocol{5});

    rmse_list = NaN(num_volts,1);
    for i = 1:num_volts
        yi = yksum(:,i);
        ymx = kcurrent_model(p, mdl_struct, protocol, volts(i), pdefault);
        yhat = sum(ymx,2);
        rmse_list(i) = sqrt(mean((yi(hold_idx+1:end) - yhat(hold_idx+1:end)).^2));
    end
    z = sum(rmse_list);
end
