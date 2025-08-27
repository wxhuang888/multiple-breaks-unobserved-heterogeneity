function beta = alpha2beta2(alpha, regime)
    p = size(alpha, 2); % number of parameters
    m = size(regime, 1) - 2; % number of regime changes
    T = regime(end) - 1; % total time periods
    N = size(alpha, 1); % number of individuals
    
    beta = zeros(N * T, p); % initialize beta matrix (N*T by p)
    
    for i = 1:N
        alphai = alpha(i,:,:); % extract alpha for individual i
        betai = zeros(T, p); % initialize beta for individual i
        
        for j = 1:m + 1
            % Assign the same alpha values to all time periods within the regime
            betai(regime(j):regime(j + 1) - 1, :) = kron(alphai(:,:, j), ones(regime(j + 1) - regime(j), 1));
        end
        
        beta((i - 1) * T + 1:i * T, :) = betai; % store beta for individual i
    end
end