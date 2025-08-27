function result = est_Group_new_emp(thisds_o, thisds, thisK,mse)
% this function get Q for c-lasso estimation

K = thisK;
unique_IDs=unique(thisds_o.N);
TT = max(thisds_o.T) - min(thisds_o.T) + 1;

% may work with TT^(-1/3)*0.05
lambda0 = mse*TT^(-1/3);
tol     = 0.05; % convergence tolerance level
R       = 80; %  maximum number of iterations



%% Initial values
p = size(thisds_o.X, 2); % Number of predictors
beta_hat0 = zeros(length(unique(thisds_o.N)), p);

for i = 1:length(unique_IDs)
    % Extract data for individual i
    yi = thisds_o.y(thisds_o.N == unique_IDs(i));
    Xi = thisds_o.X(thisds_o.N == unique_IDs(i), : );

    % Check if Xi'*Xi is invertible
    if rank(Xi' * Xi) == size(Xi, 2) % Full rank check
        beta_hat0(i, :) = (Xi' * Xi) \ (Xi' * yi); % More stable than inv()
    else
        beta_hat0(i, :) = NaN; % Mark as NaN if singular
    end
end

% Remove individuals with NaN estimates in beta_hat0
valid_ids = all(~isnan(beta_hat0), 2);
beta_hat0 = beta_hat0(valid_ids, :);
thisds_o = thisds_o(ismember(thisds_o.N, unique_IDs(valid_ids)), :);
ds_d=thisds_o;
% prepare some parameters
y = thisds_o.y;
X = thisds_o.X;
N = length(unique(ds_d.N));
%% PLS
[b_K, a, ~] = PLS_est(N, TT, y, X, beta_hat0, K, lambda0, R, tol);
[~, b, ~ , group] = report_b( b_K, a, K, p);
sum(group);
%% post estimation
NN = unique(thisds_o.N);

post_b=zeros(N,p);
est_post_std  = zeros(N,p);
est_post_std1  = zeros(N,p);
group = logical(group);
%
for i = 1:K
    this_group= group(:,i);
    g_index = NN(this_group);

    if isempty(g_index)
        est_lasso(:,i)   =  nan;
        est_post(:,i)    =  nan;
    else
        %first_none_zero = min( NN(this_group) );
        g_data = ds_d( ismember(ds_d.N, g_index), : ); % group-specific data

        %post
        Td = size(g_data,1) / size(g_index,2);
        nwestresults = nwest(g_data.y, g_data.X, floor( Td^(1/3) ));
        post_b(g_index,:)=repmat(nwestresults.beta',length(g_index),1);
        
        [vari_a_post] = var_PLS(Td, nwestresults.beta, g_data.y, g_data.X  );

        post_std= sqrt( diag( vari_a_post ) );


        est_post_std1(g_index,:)=repmat(post_std',length(g_index),1);
        
        est_post_std(g_index,:)=repmat(nwestresults.se',length(g_index),1);
        
    end
end
post_b(~any(post_b, 2), :) = [];
est_post_std(~any(est_post_std, 2), :) = [];
est_post_std1(~any(est_post_std1, 2), :) = [];
%postB = kron( post_b, ones(T,1) );
mu=zeros(N,1);
% Get unique IDs
unique_N = unique(thisds_o.N);
N = length(unique_N); % Total number of unique individuals

% Initialize output matrices
uNT_post = zeros(size(thisds.y));  % Assuming full size match
uNT_classo = zeros(size(thisds.y));

% Loop over each unique ID
for idx = 1:N
    i = unique_N(idx);  % Get the actual ID value

    % Extract y and X for individual i
    yi = thisds.y(thisds.N == i);
    Xi = thisds.X(thisds.N == i, :);
    mu(idx,1)=mean(yi-Xi*post_b(idx,:)');
    % Compute residuals after adjusting for mean
    mui_post = mean(yi - Xi * post_b(idx, :)');
    uNT_post((idx - 1) * TT + 1 : idx * TT, :) = yi - Xi * post_b(idx, :)' - mui_post;

    mui_classo = mean(yi - Xi * b(idx, :)');
    uNT_classo((idx - 1) * TT + 1 : idx * TT, :) = yi - Xi * b(idx, :)' - mui_classo;
end
classoQ =  mean(uNT_classo.^2);
post_Q =  mean(uNT_post.^2);
%% store
result.group         = group;
result.est_post = post_b;
result.est_mu = mu;
result.est_classo = b;
% % result.est_lasso     = est_lasso;
% % % result.est_lasso_std = est_lasso_std;
result.est_post_Q      = post_Q;
%result.est_post_std  = est_post_std;
result.est_post_std1 = est_post_std1;
result.est_lasso_Q   = classoQ;
% result.est_lasso_Q1   = classoQ1;
end
