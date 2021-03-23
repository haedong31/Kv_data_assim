function [bfs, bchroms] = fit1(t, y, N0, N1, N2, max_iter)
    bfs = [];
    bchroms = [];

    % initial generation
    init_gen = unifrnd(0, 0.01, N0, 28);
    init_gen(:, 9) = unifrnd(0, 1, N0, 1);
    init_gen(:, 18) = unifrnd(0, 1, N0, 1);
    init_gen(:, 27) = unifrnd(0, 1, N0, 1);
    init_gen(:, 28) = unifrnd(0, 1, N0, 1);
    cnt = 1;

    % initial evaluation
    fits = eval_fn(init_gen, t, y, N0);
    [bf, bf_idx] = min(fits);
    bchrom = init_gen(bf_idx, :);
    fprintf('Initial best fit: %f \n', bf);
%     disp(bchrom)

    bcnt = 1;
    bfs = [bfs, bf];
    bchroms = [bchroms; bchrom];

    new_gen = evolve(init_gen, fits, N0, N1, N2);

    for i=1:max_iter
        cnt = cnt + 1;
%         fprintf('\n Generation %i \n', cnt)
        
        fits = eval_fn(new_gen, t, y, N0);
        [bf, bf_idx] = min(fits);
        bchrom = new_gen(bf_idx, :);

        if (bf < bfs(bcnt))
            bcnt = bcnt + 1;
            fprintf('Best fit is updated: %f | %i \n', bf, cnt)
%           disp(bchrom)
 
            
            bfs = [bfs, bf];
            bchroms = [bchroms; bchrom];
        end
        new_gen = evolve(new_gen, fits, N0, N1, N2);
    end        
end

function [fits] = eval_fn(pop, t, y, N0)
    holdV = -70;
    [~, hold_idx] = max(y);
    P1 = 50;

    fits = zeros(1, N0);
    for i=1:N0
        A = Ktrace1(pop(i, :), holdV, hold_idx, P1, t);
        Ktrc = A(:, 1) + A(:, 2) + A(:, 3) + A(:, 4);
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
    new_gen = zeros(N0, 28);
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
