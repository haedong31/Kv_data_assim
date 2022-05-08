Surrogate-assisted calibration - fitrgp
Code argument & Model information
clc
clearvars
warning('off','all')

group = "wt";
exp_num = "test";

save_dir = strcat("calib_",exp_num,"_",group);
if ~exist(fullfile(pwd,save_dir),'dir')
    mkdir(fullfile(pwd,save_dir))
end

current_names = {'iktof','ikslow1','ikslow2','ikss'};
tune_idx = cell(5,1); % iktof, iktos, ikslow1, ikslow2, ikss
tune_idx{1} = [1, 2, 4, 5, 7, 11, 13];
tune_idx{2} = [2, 3];
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

Protocol
% voltage
holdv = -70;
minv = -30;
vstep = 10;
num_steps = 9;
ek = -91.1;

% time
ideal_holdt = 120;
ideal_endt = 4.6*1000;

% save time information later when experimental data imlported
protocol = cell(6,1);
protocol{1} = holdv;
protocol{2} = ek;
protocol{3} = minv:vstep:(minv+vstep*num_steps-1);

% meta data of experimental datasets
matching_table = readtable(fullfile(pwd,"mgat1ko_data",strcat("matching-table-",group,".xlsx")));
file_names = matching_table.trace_file_name_4half;
data_dir = fullfile("mgat1ko_data",strcat(group,"-preprocessed"));

% exclude null rows
file_names(cellfun(@isempty,file_names)) = [];
file_names = string(file_names);

Optimization loop
% for i = 1:num_files
i = 1;

% read experimental data
file_path = fullfile(pwd,data_dir,file_names(i));
trace_data = table2array(readtable(file_path));
t = trace_data(:,1);
yksum = trace_data(:,2:end);

% estiamte the critical time points
[~,ideal_hold_idx] = min(abs(t-ideal_holdt));
[~,ideal_end_idx] = min(abs(t-ideal_endt));
t = t(1:ideal_end_idx);
yksum = yksum(1:ideal_end_idx,:);

% fill out time information to complete the protocol
protocol{4} = t;
protocol{5} = t(1:ideal_hold_idx);
pulset = t(ideal_hold_idx+1:end);
protocol{6} = pulset - pulset(1);

% optimization parameters
sg = gensg(psize+3,psize);
init_dgn = scale_param(sg.design,lb,ub);
out_size = 100;
num_starts = 5;
ei_tol = sqrt(eps);

% run optimization
[sol,fval,maxei] = optim_surr(init_dgn,out_size,num_starts,ei_tol,yksum,mdl_struct,pdefault,protocol);

Custom functions
BO-routine functions
function [xnew,ynew,maxei] = optim_surr(init_dgn,out_size,num_starts,ei_tol,yksum,mdl_struct,pdefault,protocol)
    init_size = size(init_dgn);
    ts = tic;
    y = mdl_disc(init_dgn,yksum,mdl_struct,pdefault,protocol);
    te = toc(ts);

    % initial GP model
    gpmdl = fitrgp(init_dgn,y,'ComputationMethod','v',...
        'KernelFunction','ardsquaredexponential',...
        'FitMethod','exact','Optimizer','lbfgs',...
        'OptimizeHyperparameters',{'KernelFunction','KernelScale','Sigma'}, ...
        'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus',...
        'UseParallel',true,'ShowPlots',false,'Verbose',0));
    
    % sequential data acquisition
    sol = NaN(out_size,init_size(2));
    ynew = NaN(out_size,1);
    maxei = NaN(out_size-init_size(1),1);
    sol(1:init_size(1),:) = init_dgn;
    ynew(1:init_size(1),:) = y;

    for i = (init_size(1)+1):(init_size(1)+out_size)
        [xnew,ei] = ei_search(sol(1:i-1),ynew(1:i-1),gpmdl,num_starts,ei_tol);
        [running_maxei,max_idx] = max(ei);

        if isempty(max_idx)
            sol(i,:) = [];
            ynew(i) = [];
            maxei(i-init_size(1)) = [];
        else
            sol(i,:) = xnew(max_idx,:);
            ynew(i) = mdl_disc(sol(i,:),yksum,mdl_struct,pdefault,protocol);
            maxei(i-init_size(1)) = running_maxei;
            
            % update GP model
            % ts = tic;
            gpmdl  = fitrgp(sol(1:i,:),ynew(1:i),'ComputationMethod','v',...
                'KernelFunction','ardsquaredexponential',...
                'FitMethod','exact','Optimizer','lbfgs',...
                'OptimizeHyperparameters',{'KernelFunction','KernelScale','Sigma'}, ...
                'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus',...
                'UseParallel',true,'ShowPlots',false,'Verbose',0));
            % te = toc(ts);
        end
    end
end

function [xnew,ei] = ei_search(x,y,gpmdl,num_starts,ei_tol,lb,ub)
    [ymin,min_idx] = min(y);
    num_var = size(x,2);

    % multiple starting points
    if num_starts > 1
        starts = NaN(num_steps,num_var);
        starts(1,:) = x(min_idx,:);
        starts(2:end,:) = lhsdesign(num_starts-1,num_var);
    else
        starts = x(min_idx,:);
    end

    % search acquisition points by maximizing EI
    xnew = NaN(num_starts,num_var);
    ei = NaN(num_starts,1);
    obj_fun = @(xx)obj_ei(xx,ymin,gpmdl);
    optim_opts = optimoptions(@fmincon, 'Algorithm','interior-point', 'Display','off');

    for i = 1:num_starts
        if (-obj_ei(starts(i,:),ymin,gpmdl) <= ei_tol)
            ei(i) = -Inf;
        else
            [sol,fval] = fmincon(obj_fun,starts(i,:),[],[],[],[],lb,ub,[],optim_opts);
            xnew(i,:) = sol;
            ei(i) = -fval;
        end
    end

    valid_idx = ei > ei_tol;
    xnew = xnew(valid_idx,:);
    ei = ei(valid_idx);
end

function z = obj_ei(xx,fmin,gpmdl)
    [mu,s2] = predict(gpmdl,xx);
    d = fmin - mu;
    s = sqrt(s2);
    dn = d./s;
    ei = d.*normcdf(dn) + s.*normpdf(dn);
    z = -ei;
end

function d = mdl_disc(p,yksum,mdl_struct,pdefault,protocol)
    volts = protocol{3};
    hold_idx = length(protocol{5});
    yksumt = reshape(yksum,[1,size(yksum)]);
    ymx = kcurrent_model(p,mdl_struct,pdefault,protocol,volts);
    yksum_hat = squeeze(sum(ymx,4));
    rmse_list = squeeze(sqrt(mean((yksum_hat(:,hold_idx+1:end,:)-yksumt(:,hold_idx+1:end,:)).^2,2)));
    d = sum(rmse_list,2);
end
Helper functions
function [mdl_struct, psize] = gen_mdl_struct(current_names, tune_idx)
    num_currents = length(current_names);
    
    % index 1
    idx_info1 = cell(1,num_currents);
    for i = 1:num_currents
        switch current_names{i}
            case 'iktof'
                idx_info1{i} = tune_idx{1};
            case 'iktos'
                idx_info1{i} = tune_idx{2};
            case 'ikslow1'
                idx_info1{i} = tune_idx{3};
            case 'ikslow2'
                idx_info1{i} = tune_idx{4};
            case 'ikss'
                idx_info1{i} = tune_idx{5};
        end
    end

    % index 2
    idx_info2 = cell(1,num_currents);
    psize = 0;
    for i = 1: num_currents
        psize_old = psize;
        psize = psize + length(idx_info1{i});
        idx_info2{i} = (1+psize_old):(psize);
    end

    model_info = [current_names; idx_info1; idx_info2];
    field_names = ["name","idx1","idx2"];
    mdl_struct = cell2struct(model_info, field_names, 1);    
end

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
        scaledp(i, :) = unitp(i,:).*(ub-lb)' + lb';
    end    
end
Current models
function ymx = kcurrent_model(p,mdl_struct,pdefault,protocol,volts)
    num_pts = size(p,1);
    current_names = {mdl_struct.name};
    idx_info1 = {mdl_struct.idx1};
    idx_info2 = {mdl_struct.idx2};
    
    % define param_kslow1 first just in case it would be used in advance
    matching_idx = strcmp(current_names, 'ikslow1');
    param_kslow1 = repmat(pdefault{3},num_pts,1);
    if any(matching_idx)
        param_kslow1(:,idx_info1{matching_idx}) = p(:,idx_info2{matching_idx});
    end
    
    ymx = zeros(num_pts, length(protocol{4}), length(volts),length(current_names));
    for i = 1:length(current_names)
        switch current_names{i}
            case 'iktof'
                param = repmat(pdefault{1},num_pts,1);
                param(:,idx_info1{i}) = p(:,idx_info2{i});
                ymx(:,:,:,i) = iktof(param, protocol, volts);
            case 'iktos'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);
            
                param = NaN(num_pts,num_param);
                param(:,shared_idx) = param_kslow1(:,shared_idx);
                uniq_param = repmat(pdefault{2},num_pts,1);
                uniq_param(:,idx_info1{i}) = p(:,idx_info2{i});
                param(:,uniq_idx) = uniq_param;
                ymx(:,:,i) = iktos(param, protocol, volts);
            case 'ikslow1'
                ymx(:,:,:,i) = ikslow1(param_kslow1, protocol, volts);
            case 'ikslow2'
                num_param = 11;
                shared_idx = [1:7,9];
                uniq_idx = setdiff(1:num_param, shared_idx);
            
                param = NaN(num_pts,num_param);
                param(:,shared_idx) = param_kslow1(:,shared_idx);
                uniq_param = repmat(pdefault{4},num_pts,1);
                uniq_param(:,idx_info1{i}) = p(:,idx_info2{i});
                param(:,uniq_idx) = uniq_param;
                ymx(:,:,:,i) = ikslow2(param, protocol, volts);
            case 'ikss'
                num_param = 7;
                shared_idx1 = 1:3;
                shared_idx2 = [1,3,4];
                uniq_idx = setdiff(1:num_param, shared_idx1);
            
                param = NaN(num_pts,num_param);
                param(:,shared_idx1) = param_kslow1(:,shared_idx2);
                uniq_param = repmat(pdefault{5},num_pts,1);
                uniq_param(:,idx_info1{i}) = p(:,idx_info2{i});
                param(:,uniq_idx) = uniq_param;
                ymx(:,:,:,i) = ikss(param, protocol, volts);
        end
    end
end

function y = iktof(p, protocol, volts)
    gmax = p(:,13);
    act0 = 0.4139033547e-02;
    inact0 = 0.9999623535e+00;
    
    holdv = protocol{1};
    holdv = ones(1,length(volts))*holdv;
    num_volts = length(volts);
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};

    y = zeros(size(p,1),length(t),num_volts);

    % current equation at holding volt
    df = reshape(holdv-ek,[1,1,num_volts]);
    kv_hold = ikto_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(:,:,1), kv_hold(:,:,3));
    inact_hold = hh_model(holdt, inact0, kv_hold(:,:,2), kv_hold(:,:,4));
    y(:,1:hold_idx,:) = gmax.*(act_hold.^3).*(inact_hold).*df;

    % current equation at pulse volt
    df = reshape(volts-ek,[1,1,num_volts]);
    kv_pulse = ikto_kinetic_variables(p, volts);
    act_pulse = hh_model(pulset, act0, kv_pulse(:,:,1), kv_pulse(:,:,3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(:,:,2), kv_pulse(:,:,4));
    y(:,(hold_idx+1):end,:) = gmax.*(act_pulse.^3).*(inact_pulse).*df;
end

function y = iktos(p, protocol, volts)
    gmax = p(:,11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    holdv = protocol{1};
    holdv = ones(1,length(volts))*holdv;
    num_volts = length(volts);
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(size(p,1),length(t),num_volts);

    % current equation at holding
    df = reshape(holdv-ek,[1,1,num_volts]);
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(:,:,1), kv_hold(:,:,3));
    inact_hold = hh_model(holdt, inact0, kv_hold(:,:,2), kv_hold(:,:,4));
    y(:,1:hold_idx,:) = gmax.*(act_hold.*inact_hold).*df;
    
    % current equation at pulse voltage
    df = reshape(volts-ek,[1,1,num_volts]);
    kv_pulse = rectifier_kinetic_variables(p, volts);
    act_pulse = hh_model(pulset, act0, kv_pulse(:,:,1), kv_pulse(:,:,3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(:,:,2), kv_pulse(:,:,4));
    y(:,(hold_idx+1):end,:) = gmax.*(act_pulse.*inact_pulse).*df;
end

function y = ikslow1(p, protocol, volts)
    gmax = p(:,11);
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;
    
    holdv = protocol{1};
    holdv = ones(1,length(volts))*holdv;
    num_volts = length(volts);
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(size(p,1),length(t),num_volts);

    % current equation at holding
    df = reshape(holdv-ek,[1,1,num_volts]);
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(:,:,1), kv_hold(:,:,3));
    inact_hold = hh_model(holdt, inact0, kv_hold(:,:,2), kv_hold(:,:,4));
    y(:,1:hold_idx,:) = gmax.*(act_hold.*inact_hold).*df;
    
    % current equation at pulse voltage
    df = reshape(volts-ek,[1,1,num_volts]);
    kv_pulse = rectifier_kinetic_variables(p, volts);
    act_pulse = hh_model(pulset, act0, kv_pulse(:,:,1), kv_pulse(:,:,3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(:,:,2), kv_pulse(:,:,4));
    y(:,(hold_idx+1):end,:) = gmax.*(act_pulse.*inact_pulse).*df;
end

function y = ikslow2(p, protocol, volts)
    gmax = p(:,11); 
    act0 = 0.5091689794e-03;
    inact0 = 0.9980927689;

    holdv = protocol{1};
    holdv = ones(1,length(volts))*holdv;
    num_volts = length(volts);
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(size(p,1),length(t),num_volts);

    % current equation at holding
    df = reshape(holdv-ek,[1,1,num_volts]);
    kv_hold = rectifier_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(:,:,1), kv_hold(:,:,3));
    inact_hold = hh_model(holdt, inact0, kv_hold(:,:,2), kv_hold(:,:,4));
    y(:,1:hold_idx,:) = gmax.*(act_hold.*inact_hold).*df;
    
    % current equation at pulse voltage
    df = reshape(volts-ek,[1,1,num_volts]);
    kv_pulse = rectifier_kinetic_variables(p, volts);
    act_pulse = hh_model(pulset, act0, kv_pulse(:,:,1), kv_pulse(:,:,3));
    inact_pulse = hh_model(pulset, inact0, kv_pulse(:,:,2), kv_pulse(:,:,4));
    y(:,(hold_idx+1):end,:) = gmax.*(act_pulse.*inact_pulse).*df;
end

function y = ikss(p, protocol, volts)
    gmax = p(:,7);
    act0 = 0.5091689794e-03;

    holdv = protocol{1};
    holdv = ones(1,length(volts))*holdv;
    num_volts = length(volts);
    ek = protocol{2};
    t = protocol{4};
    holdt = protocol{5};
    hold_idx = length(holdt);
    pulset = protocol{6};
    
    y = zeros(size(p,1),length(t),num_volts);

    % current equation at holding
    df = reshape(holdv-ek,[1,1,num_volts]);
    kv_hold = ikss_kinetic_variables(p, holdv);
    act_hold = hh_model(holdt, act0, kv_hold(:,:,1), kv_hold(:,:,2));
    y(:,1:hold_idx,:) = gmax.*act_hold.*df;

    % current equation at pulse voltage
    df = reshape(volts-ek,[1,1,num_volts]);
    kv_pulse = ikss_kinetic_variables(p, volts);
    act_pulse = hh_model(pulset, act0, kv_pulse(:,:,1), kv_pulse(:,:,2));
    y(:,(hold_idx+1):end,:) = gmax.*act_pulse.*df;
end

function kv = ikto_kinetic_variables(p,volts)
    kv = NaN(size(p,1),length(volts),8);
    
    alpha1 = p(:,7).*exp(p(:,5).*(volts+p(:,1)));
    beta1 = p(:,8).*exp(-p(:,6).*(volts+p(:,1)));
    
    alpha2_temp1 = p(:,9).*exp((volts+p(:,2))./(-1.0.*p(:,4)));
    alpha2_temp2 = p(:,10).*exp((volts+p(:,2)+p(:,3))./(-1.0.*p(:,4)));
    alpha2 = alpha2_temp1./(1.0+alpha2_temp2);

    beta2_temp1 = p(:,11).*exp((volts+p(:,2)+p(:,3))./p(:,4));
    beta2_temp2 = p(:,12).*exp((volts+p(:,2)+p(:,3))./p(:,4));
    beta2 = beta2_temp1./(1.0+beta2_temp2);

    kv(:,:,1) = alpha1./(alpha1+beta1);
    kv(:,:,2) = alpha2./(alpha2+beta2);
    kv(:,:,3) = 1./(alpha1+beta1);
    kv(:,:,4) = 1./(alpha2+beta2);
    kv(:,:,5) = alpha1;
    kv(:,:,6) = beta1;
    kv(:,:,7) = alpha2;
    kv(:,:,8) = beta2;
end

function kv = rectifier_kinetic_variables(p,volts)
    kv = NaN(size(p,1),length(volts),4);
    kv(:,:,1) = 1.0./(1.0+exp(-(p(:,1)+volts)./p(:,4))); % ass
    kv(:,:,2) = 1.0./(1.0+exp((p(:,2)+volts)./p(:,5))); % iss
    kv(:,:,3) = p(:,7)./(exp(p(:,6).*(volts+p(:,3))) + exp(-p(:,6).*(volts+p(:,3))))+p(:,9); % taua
    kv(:,:,4) = p(:,10) - p(:,8)./(1.0+exp((p(:,2)+volts)./p(:,5))); % taui
end

function kv = ikss_kinetic_variables(p,volts)
    kv = NaN(size(p,1),length(volts),2);
    kv(:,:,1) = 1.0./(1.0+exp(-(p(:,1)+volts)./p(:,3)));
    kv(:,:,2) = p(:,5)./(exp(p(:,4).*(volts+p(:,2))) + exp(-p(:,4).*(volts+p(:,2)))) + p(:,6);
end

function y = hh_model(t, ss0, ss, tau)
    tau = reshape(tau,[size(tau,1),1,size(tau,2)]);
    ss = reshape(ss,[size(ss,1),1,size(ss,2)]);
    y = ss - (ss-ss0).*exp(-t'./tau);
end
