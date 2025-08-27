function alpha = simu_genalpha_N1F_RS(sigma_a, jmin, jmax, p, N, m0, N1B)
    % Step 1: Initialize alpha0 as an (N, p, m0+1) matrix
    alpha0 = zeros(N, p, m0 + 1); 
    
    % Step 2: Generate alpha0 for each individual
    for i = 1:N
        % First regime: Initialize the first slice of alpha0 for individual i
        alpha0(i, :, 1) = sigma_a * randn(1, p);
        
        % Generate subsequent regimes for individual i
        for j = 2:m0 + 1
            alpha0(i, :, j) = alpha0(i, :, j - 1) + sign(randn(1, p)) .* (jmin + (jmax - jmin) * rand(1, p));
        end
    end
    
    % Step 3: Initialize alpha
    alpha = zeros(N, p, m0 + 1); % Initialize alpha (N x p x (m0+1))
    
    % Step 4: Randomly select N1B individuals to remain unchanged across all regimes
    fixed_individuals = randperm(N, N1B); % Randomly select N1B indices
    changing_individuals = setdiff(1:N, fixed_individuals); % Remaining individuals will change
    
    % Assign values to alpha
    % For selected individuals in fixed_individuals, all regimes take the same value
    for idx = fixed_individuals
        alpha(idx, :, :) = repmat(alpha0(idx, :, 1), [1, 1, m0 + 1]);
    end
    
    % For remaining individuals in changing_individuals, allow changes across regimes
    for idx = changing_individuals
        alpha(idx, :, :) = alpha0(idx, :, :);
    end
end