function V = getV(ds, thiss1, thiss2, thisT2, thisdT, thislambda0)

thisT = thiss2 - thiss1 + 1;
unique_N = unique(ds.N);  % Get unique values of ds.N
thisN = size(unique_N, 1);  % Number of unique groups
thisp = size(ds.X, 2);  % Dimension of parameters

V = nan(thiss2, 1);  % Initialize V as NaN
for thistt = thiss1 + thisT2 + thisp : thisdT : thiss2 - thisT2 - thisp  % Searching range
    % thistt
    res_pre = zeros(thisN,1);
    res_post = zeros(thisN,1);
    for i = 1:thisN
        current_N = unique_N(i);  % Get the current value of ds.N being processed
        
        % Data
        g_data_pre_i = ds(ds.N == current_N & ds.T >= thiss1 & ds.T < thistt, :);
        g_data_post_i = ds(ds.N == current_N & ds.T >= thistt & ds.T <= thiss2, :);

        % Add a constant (intercept) to the design matrices X
        %g_data_pre_i.X = [ones(size(g_data_pre_i.X, 1), 1), g_data_pre_i.X];
        %g_data_post_i.X = [ones(size(g_data_post_i.X, 1), 1), g_data_post_i.X];

        % Check if g_data_pre_i.y or g_data_post_i.y is a zero vector
        if all(g_data_pre_i.y == 0)
            res_pre(i, :) = 0;  % If g_data_pre_i.y is a zero vector, set res_pre to an empty array
        else
            % Estimate \tilde{\beta}^{a} and calculate residuals
            beta_prei = regress(g_data_pre_i.y, g_data_pre_i.X);
            res_prei = g_data_pre_i.y - g_data_pre_i.X * beta_prei;
            res_pre(i, :)=sum(res_prei.^2);
        end
        
        if all(g_data_post_i.y == 0)
            res_post(i, :) = 0;  % If g_data_post_i.y is a zero vector, set res_post to an empty array
        else
            % Estimate \tilde{\beta}^{b} and calculate residuals
            beta_posti = regress(g_data_post_i.y, g_data_post_i.X);
            res_posti= g_data_post_i.y - g_data_post_i.X * beta_posti;
            res_post(i, :)=sum(res_posti.^2);
        end        
    end

    Q_pre = sum(res_pre);
    Q_post = sum(res_post);

    % Calculate V
    V(thistt, 1) = Q_pre + Q_post;
    if thistt ~= thiss1 && thistt ~= thiss2
        V(thistt, 1) = V(thistt, 1) + thislambda0;
    end
    
end

end
