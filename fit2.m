function [bfs, bchroms] = fit2(time_space, y, N0, N1, N2, max_iter)
    bfs = [];
    bchroms = [];

    % time_space
    t = time_space{1};
    tH = time_space{2};
    tP1_adj = time_space{3};

    % initial population; Ito
    init_Ito = zeros(N0, 12);
    init_Ito(:, 1:3) = unifrnd(0, 70, N0, 3);
    init_Ito(:, 4) = unifrnd(0.1, 15, N0, 1);
    init_Ito(:, 5:12) = unifrnd(0.0001, 0.5, N0, 8);

    % initial population; Ito
    init_IKslow1 = zeros(N0, 11);
    init_IKslow1(:, 1) = unifrnd(0, 50, N0, 1);
    init_IKslow1(:, 2) = unifrnd(0.1, 15, N0, 1);
    init_IKslow1(:, 3:5) = unifrnd(0.1, 5, N0, 3);
    init_IKslow1(:, 6) = unifrnd(10, 80, N0, 1);
    init_IKslow1(:, 7) = unifrnd(0.1, 15, N0, 1);
    init_IKslow1(:, 8) = unifrnd(1500, 5000, N0, 1);
    init_IKslow1(:, 9) = unifrnd(100, 1000, N0, 1);
    init_IKslow1(:, 10) = unifrnd(0, 50, N0, 1);
    init_IKslow1(:, 11) = unifrnd(0.1, 15, N0, 1);

    % initial population; Ito
    init_IKslow2 = zeros(N0, 11);
    init_IKslow2(:, 1) = unifrnd(0, 50, N0, 1);
    init_IKslow2(:, 2) = unifrnd(0.1, 15, N0, 1);
    init_IKslow2(:, 3:5) = unifrnd(0.1, 5, N0, 3);
    init_IKslow2(:, 6) = unifrnd(10, 80, N0, 1);
    init_IKslow2(:, 7) = unifrnd(0.1, 15, N0, 1);
    init_IKslow2(:, 8) = unifrnd(5000, 15000, N0, 1);
    init_IKslow2(:, 9) = unifrnd(100, 1000, N0, 1);
    init_IKslow2(:, 10) = unifrnd(0, 50, N0, 1);
    init_IKslow2(:, 11) = unifrnd(0.1, 15, N0, 1);

    % initial population; Ito
    init_Iss = unifrnd(0.1, 5, N0, 1);

    % initial population; Ito
    initG = zeros(N0, 3);
    initG(:, 1:3) = unifrnd(0.1, 1, N0, 3);

    init_pop = [init_Ito, init_IKslow1, init_IKslow2, init_Iss, initG];
    cnt = 1;
    
    fits = eval_fn(init_pop, y, time_space, N0);
    [bf, bf_idx] = min(fits);
    bchrom = init_pop(bf_idx, :);
    fprintf('Initial best fit: %f \n', bf);

    bcnt = 1;
    bfs = [bfs, bf];
    bchroms = [bchroms; bchrom];

    new_gen = evolve(init_pop, fits, N0, N1, N2);

    for i=1:max_iter
        cnt = cnt + 1;
        % fprintf('\n Generation %i \n', cnt)
        
        fits = eval_fn(new_gen, y, time_space, N0);
        [bf, bf_idx] = min(fits);
        bchrom = new_gen(bf_idx, :);

        if (bf < bfs(bcnt))
            bcnt = bcnt + 1;
            
            if (mod(bcnt, 10) == 0)
                fprintf('Best fit is updated: %f | %i \n', bf, cnt)
                % disp(bchrom)
            end
            
            bfs = [bfs, bf];
            bchroms = [bchroms; bchrom];
        end
        new_gen = evolve(new_gen, fits, N0, N1, N2);
    end            
end

function [fits] = eval_fn(pop, y, time_space, N0)
    t = time_space{1};
    fits = zeros(1, N0);
    for i=1:N0
        param = cell(1, 5);
        param{1} = pop(i, 1:12);
        param{2} = pop(i, 13:23);
        param{3} = pop(i, 24:34);
        param{4} = pop(i, 35);
        param{5} = pop(i, 36:38);
        k = Ktrace2(param, -70, 50, time_space);
        Ktrc = k(:, 1) + k(:, 2) + k(:, 3) + k(:, 4);
        fits(i) = sqrt((1/length(t))*sum((y-Ktrc).^2));
    end
end

function [new_gen] = evolve(pop, fits, N0, N1, N2)
    [~, elite_idx] = mink(fits, N1);
    elites = pop(elite_idx, :);

    mean_elite = mean(elites, 1);
%     sigs = std(elites);
    elites(end, :) = mean_elite;

    cnt = 1;
    new_gen = zeros(N0, 38);
    new_gen(1:N1, :) = elites;
    for i=1:N1
        elite = elites(i, :);
        for j=1:N2
            offspr = elite + normrnd(0, 0.01);
            offspr = abs(offspr);
            new_gen(N1+cnt, :) = offspr;
            cnt = cnt + 1;
        end
    end
end
