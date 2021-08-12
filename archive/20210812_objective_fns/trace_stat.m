function [stats] = trace_stat(t, current_trace)
    % statistics of current trace
    % stat 1: peak time
	% stat 2: peak current value
	% stat 3: tau, peak reduced by exp(-1) (~63%)
	% stat 4: current at 1/3 of interval [peak, tau]
	% stat 5: current at 2/3 of interval [peak, tau]
	% stat 6: current at last
    
    stats = zeros(6, 1);

    % truncate trace
    [peak, peak_idx] = max(current_trace);
    current_trace_trunc = current_trace(peak_idx:end);
    t_trunc = t(peak_idx:end);
    t_trunc = t_trunc - t_trunc(1);
    
    % stat 1
    stats(1) = t(peak_idx);
    
    % stat 2
    stats(2) = peak;
    
    % stat 3
    [~, tau_idx] = min(abs(peak*exp(-1) - current_trace_trunc));
    stats(3) = t_trunc(tau_idx);

    % stat 4
    [~, stat45_jump_size] = min(abs(t_trunc(tau_idx)/3 - t_trunc));
    stats(4) = current_trace_trunc(stat45_jump_size);

    % stat 5
    stats(5) = current_trace_trunc(stat45_jump_size*2);

    % stat 6
    stats(6) = current_trace(end);
end
