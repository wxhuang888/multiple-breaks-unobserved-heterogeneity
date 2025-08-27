function result = est_Group2_lambda(thisds, ds, thisK, lambda0)
% this function get Q for c-lasso estimation

K = thisK;
N = max(thisds.N);
T = max(thisds.T) - min(thisds.T) + 1;

%lambda0 = T^(-1/4);
tol     = 0.0001; % convergence tolerance level
R       = 80; %  maximum number of iterations        
        
% prepare some parameters
y = thisds.y;
X = thisds.X;
p = size(X,2);

%% initial values
beta_hat0 = zeros(N, p);
for i = 1:N
    yi = thisds.y(thisds.N == i );
    Xi = thisds.X(thisds.N == i, : );
    beta_hat0(i,:) = ( inv(Xi'*Xi)*Xi'*yi )';
end

%% PLS
[b_K, a, ~] = PLS_est(N, T, y, X, beta_hat0, K, lambda0, R, tol);
[~, b, ~ , group] = report_b( b_K, a, K, p);
%classoB = kron( b, ones(T,1) );
for i = 1:N
    yi = ds.y(ds.N == i);
    Xi = ds.X(ds.N == i, : );
       
    % demean
    yid = yi-mean(yi);
    Xid = Xi-mean(Xi);
    
    y(ds.N==i) = yid;
    X(ds.N == i, :) =Xid;
end
% 
code=ds.N;
period=ds.T;
ds_d = dataset( code, period, y, X );
ds_d.Properties.VarNames = {'N'  'T'  'y'  'X'};
%% post estimation
NN = 1:N;

est_lasso     = zeros(K,p);
est_lasso_std = zeros(K,p);
est_post      = zeros(K,p);
est_post_std  = zeros(K,p);
% post_b=zeros(N,p);
group = logical(group);

% 
for i = 1:K
    this_group= group(:,i);
    g_index = NN(this_group);

    if isempty(g_index)
        est_lasso(g_index,:)   =  nan;
        est_lasso_std(g_index,:)   =  nan;
        est_post(g_index,:)    =  nan;
        est_post_std(g_index,:)   =  nan;
    else
        first_none_zero = min( NN(this_group) );
        g_data_o = thisds( ismember(thisds.N, g_index), : ); % group-specific data
        classo_a = b(first_none_zero,:);
        est_lasso(i,:)=classo_a;
       % [~,classo] = post_est( N, classo_a, g_data, this_group);
        [vari_a] = var_PLS(T, classo_a',  g_data_o.y,  g_data_o.X  );
        classo.std=sqrt(vari_a);
        est_lasso_std(i,:)= classo.std;
        
        g_data = ds_d( ismember(ds_d.N, g_index), : );

        %post-nw
        nwestresults = nwest(g_data.y, g_data.X, floor( 1.3*T^(1/2) ) );   
        est_post(i,:)= nwestresults.beta;
        est_post_std(i,:)=nwestresults.se;

    end

end

%% store
result.group         = group;
result.est_lasso     = est_lasso;
result.est_lasso_std = est_lasso_std;
result.est_post      = est_post;
result.est_post_std  = est_post_std;

end
