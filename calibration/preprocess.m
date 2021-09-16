%% data preprocessing: remove upward spikes in preprocess2
clc
close all
clear variables

group = 'ko';

matching_table = readtable(fullfile(pwd, 'data', strcat('matching-table-', group, '.xlsx')));
file_names = matching_table.trace_file_name_4half;

% exclude row not having 4.5-sec data
loop_idx = [];
[num_files, ~] = size(matching_table);
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

% manually trim data inside for loop
for l = 1:length(loop_idx)
    i = loop_idx(l);
    fprintf('[File %i/%i] %s \n', l, length(loop_idx), file_names{i})
    read_path = fullfile(pwd, 'data', strcat(group, '-preprocessed2'), file_names{i});
    save_path = fullfile(pwd, 'data', strcat(group, '-preprocessed'), file_names{i});
    
    trace_data = table2array(readtable(read_path));
    yksum = trace_data(:, 2:end);
    [~, num_volts] = size(yksum);
    
    yksum_copy = yksum;
    for j = 1:num_volts
        y = yksum_copy(:, j);
        y(y <0 ) = 0;
        yksum_copy(:, j) = y;
    end
    trace_data(:, 2:end) = yksum_copy;
    writematrix(trace_data, save_path)
end

%% data preprocessing: Mgat1KO remove noise from start and end of current
clc
close all
clear variables

matching_table = readtable('./data/matching-table-ko.xlsx');
file_names = matching_table.trace_file_name_4half;

% exclude row not having 4.5-sec data
loop_idx = [];
[num_files, ~] = size(matching_table);
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

% manually trim data inside for loop
for i = loop_idx
    read_path = fullfile(pwd, 'data', 'ko-preprocessed', file_names{i});
    save_path = fullfile(pwd, 'data', 'ko-preprocessed2', file_names{i});
    
    trace_data = table2array(readtable(read_path));
    t = trace_data(:, 1);
    yksum = trace_data(:, 2:end);
    [~, num_volts] = size(yksum);
    
    figure(1)
    plot(t, yksum(:, 1))
    hold on
    yline(0, 'Color','red', 'LineWidth',2);
    for j = 2:num_volts
        plot(t, yksum(:, j))
    end
    hold off
    
    yksum_copy = yksum;
    for j = 1:num_volts
        y = yksum_copy(:, j);
        [~, early_idx] = min(abs(t - 130));
        yearly = y(1:early_idx);
        
        yearly(yearly > 0.6) = 0;
        y(1:early_idx) = yearly;
        yksum_copy(:,j) = y;
        
        figure(j+1)
        plot(t, y)
    end
        
    trace_data(:, 2:end) = yksum_copy;
    writematrix(trace_data, save_path);
end

%% data preprocessing: WT remove noise from start and end of current
clc
close all
clear variables

matching_table = readtable('./data/matching-table-wt.xlsx');
file_names = matching_table.trace_file_name_4half;

% exclude row not having 4.5-sec data
loop_idx = [];
[num_files, ~] = size(matching_table);
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

% manually trim data inside for loop
for i = loop_idx
    read_path = fullfile(pwd, 'data', 'wt-preprocessed', file_names{i});
    save_path = fullfile(pwd, 'data', 'wt-preprocessed2', file_names{i});
    
    trace_data = table2array(readtable(read_path));
    t = trace_data(:, 1);
    yksum = trace_data(:, 2:end);
    [~, num_volts] = size(yksum);
    
    figure(1)
    plot(t, yksum(:, 1))
    hold on
    yline(0, 'Color','red', 'LineWidth',2);
    for j = 2:num_volts
        plot(t, yksum(:, j))
    end
    hold off
    
    yksum_copy = yksum;
    for j = 1:num_volts
        y = yksum_copy(:, j);
        [~, early_idx] = min(abs(t - 130));
        yearly = y(1:early_idx);
        
        yearly(yearly > 0.6) = 0;
        y(1:early_idx) = yearly;
        yksum_copy(:,j) = y;
        
        figure(j+1)
        plot(t, y)
    end
        
    trace_data(:, 2:end) = yksum_copy;
    writematrix(trace_data, save_path);
end

%% data preprocessing: downsample & normalization
clc
close all
clear variables

matching_table = readtable('./data/matching-table-wt.xlsx');
file_names = matching_table.trace_file_name_4half;
caps = matching_table.cap;

% exclude row not having 4.5-sec data
loop_idx = [];
[num_files, ~] = size(matching_table);
for i = 1:num_files
    if isempty(file_names{i})
        continue
    end
    loop_idx = [loop_idx, i];
end

ideal_end_time = 4.6*1000;
for i = loop_idx
    read_path = fullfile(pwd, 'data', 'wt-preprocessed', file_names{i});
    save_path = fullfile(pwd, 'data', 'wt-preprocessed', file_names{i});

    % read data
    trace_data = table2array(readtable(read_path));

    % downsample
    trace_data = downsample(trace_data, 20);

    % normalize
    trace_data(:, 2:end) = trace_data(:, 2:end) ./ caps(i); 

    % estimate time points
    [~, ideal_end_idx] = min(abs(trace_data(:, 1) - ideal_end_time));
    
    % cut data
    trace_data = trace_data(1:ideal_end_idx, :);

    % save data
    writematrix(trace_data, save_path);
end

%% preprocessing for wt-4.5 / 12.10 / 10/29/2015 Cell 3 / 15o29009
clc
close all
clear variables

data = readtable('./data/wt/15o29009.xlsx');
col_names = data.Properties.VariableNames;
data = table2array(data);

% downsample
data_sampled = downsample(data, 20);
t = data_sampled(:, 1);
[~, num_volts] = size(data_sampled);
num_volts = num_volts - 1;

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

% remove negative values
for i = 1:num_volts
    y = data_sampled(:,i+1);
    y(y < 0) = 0;
    data_sampled(:,i+1) = y;
end

% chop off tail
ideal_hold_time = 120;
[~, hold_idx] = min(abs(t - ideal_hold_time));
ideal_end_time = ideal_hold_time + 4500;
[~, tail_idx] = min(abs(t - ideal_end_time));
data_sampled = data_sampled(1:tail_idx, :);
t = data_sampled(:, 1);

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

% each current
volt_idx = 11;
y = data_sampled(:, 1+volt_idx);
plot(t, y)

impute_val = mean(y(1:hold_idx));
y(y > 16000) = impute_val;
plot(t, y)

data_sampled(:, 1+volt_idx) = y;

% normalization
cap = 227.00;
data_sampled(:, 2:end) = data_sampled(:, 2:end)./cap;

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

data_processed = array2table(data_sampled, 'VariableNames',col_names);
writetable(data_processed, './data/wt-preprocessed/15o29009.xlsx')

%% preprocessing for wt-25 / 12.10 / 10/29/2015 Cell 3 / 15o29010


%% preprocessing for ko-4.5 / 15.40 / 11/24/2015 Cell 5 / 15n24005.xlsx
clc
close all
clear variables

data = readtable('./data/ko/15n24005.xlsx');
col_names = data.Properties.VariableNames;
data = table2array(data);

% downsample
data_sampled = downsample(data, 20);
t = data_sampled(:, 1);
[~, num_volts] = size(data_sampled);
num_volts = num_volts - 1;

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

% remove negative values
for i = 1:num_volts
    y = data_sampled(:,i+1);
    y(y < 0) = 0;
    data_sampled(:,i+1) = y;
end

% chop off tail
ideal_hold_time = 120;
[~, hold_idx] = min(abs(t - ideal_hold_time));
ideal_end_time = ideal_hold_time + 4500;
[~, tail_idx] = min(abs(t - ideal_end_time));
data_sampled = data_sampled(1:tail_idx, :);
t = data_sampled(:, 1);

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

% each current
volt_idx = 11;
y = data_sampled(:, 1+volt_idx);
plot(t, y)

impute_val = mean(y(1:hold_idx));
y(y > 5000) = impute_val;
plot(t, y)

data_sampled(:, 1+volt_idx) = y;

% normalization
cap = 177.00;
data_sampled(:, 2:end) = data_sampled(:, 2:end)./cap;

% visualization
plot(t, data_sampled(:,1+1))
hold on
for i=2:num_volts
    plot(t, data_sampled(:,i+1))
end
hold off

data_processed = array2table(data_sampled, 'VariableNames',col_names);
writetable(data_processed, './data/ko-preprocessed/15n24005.xlsx')

%% preprocessing for ko-25 / 15.40 / 11/24/2015 Cell 5 / 15n24014


%% preprocessing for 4.5-2-avg-ko data
clc
close all
clear variables

exp_ksum = readtable('./4.5s-avg-ko-orig.csv');
col_names = exp_ksum.Properties.VariableNames;

exp_ksum = table2array(exp_ksum);

t = exp_ksum(:,1);
[~, num_volts] = size(exp_ksum);
num_volts = num_volts - 1;

% remove sharp negative peaks at the early phase
for i=1:num_volts
    y = exp_ksum(:,i+1);
    y(y < 0) = 0;
    exp_ksum(:,i+1) = y;
end

% visualize experimental data 
plot(t, exp_ksum(:,1+1))
hold on
for i=2:num_volts
    plot(t, exp_ksum(:,i+1))
end
hold off

% remove sharp pulse in -50 ~ 0 mV
% -50 mV
y = exp_ksum(:,6);
[peak, peak_idx] = max(y);
impute_val = mean(y(1:100));
y(y>0.5) = impute_val;
y(255:257) = impute_val;
plot(t, y)
xlabel('Time (ms)')
ylabel('Current (pA/pF)')
exp_ksum(:,2) = y;

% - 40 mV;
y = exp_ksum(:,3);
[peak, peak_idx] = max(y);
impute_val = mean(y(1:155));
y(y>0.5) = impute_val;
y(250:255) = impute_val;
plot(t,y)
exp_ksum(:,3) = y;

% -30 mV;
y = exp_ksum(:,4);
[peak, peak_idx] = max(y);
impute_val = mean(y(1:180));
y(y>0.6) = impute_val;
plot(t,y)
exp_ksum(:,4) = y;

% -20 mV
y = exp_ksum(:,5);
[peak, peak_idx] = max(y);
impute_val = mean(y(1:200));
y(y>1.5) = impute_val;
plot(t,y)
exp_ksum(:,5) = y;

% -10 mV
y = exp_ksum(:,6);
[peak, peak_idx] = max(y);
impute_val = mean(y(1:220));
y(y>3.5) = impute_val;
plot(t,y)
exp_ksum(:,6) = y;

% save
exp_ksum = array2table(exp_ksum, 'VariableNames',col_names);
writetable(exp_ksum,'./4.5s-avg-ko.csv');

%% preprocessing for 4.5-2-avg-wt data
clc
close all
clear variables

ds = readtable('./4.5s-avg-wt-orig.csv');
col_names = ds.Properties.VariableNames;

ds = table2array(ds);

t = ds(:,1);
traces = ds(:,2:end);
[~, num_volts] = size(traces);

% remove negative pulses in early phase
for i=1:num_volts
    y = traces(:,i);
    y(y < 0) = 0;
    traces(:,i) = y;
end

% cut tails
tail_idx = find(t == 4620);
t(tail_idx:end) = [];
traces(tail_idx:end,:) = [];

% remove sharp pulse in -50 ~ 0 mV
% -50 mV
y = traces(:,1);
plot(t, y)
impute_val = mean(y(1:300));
y(y>0.4) = impute_val;
y(310:315) = impute_val;
traces(:,1) = y;

% -40 mV
y = traces(:,2);
plot(t,y)
impute_val = mean(y(1:300));
y(y>0.4) = impute_val;
traces(:,2) = y;

% -30 mV
y = traces(:,3);
plot(t, y)
impute_val = mean(y(1:300));
y(y>1.5) = impute_val;
y(314:316) = impute_val;
traces(:,3) = y;

% visualize all traces
plot(t, traces(:,1))
hold on
for i=1:num_volts
    plot(t, traces(:,i))
end
hold off
xlabel('Time (ms)')
ylabel('Current (pA/pF)')

new_ds = [t, traces];
new_ds = array2table(new_ds, 'VariableNames',col_names);
writetable(new_ds,'./4.5s-avg-wt.csv')
