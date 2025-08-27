function result = est_Group_new(thisds, ds, thisK)
% this function get Q for c-lasso estimation

K = thisK;
N = max(thisds.N);
T = max(thisds.T) - min(thisds.T) + 1;

lambda0 = T^(-1/4);
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


% uNT=zeros(N*T,1);
% for i=1:N
% yi=ds.y(ds.N == i );
% Xi=ds.X(ds.N == i );
% % yi(1,:)=[];
% % Xi(1,:)=[];
% mu=mean(yi-Xi*b(i,:)');
% uNT((i-1)*T+1:i*T,:)=yi-Xi*b(i,:)'-mu;
% end    
% classoQ1 =  mean(uNT.^2);

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

post_b=zeros(N,p);
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
        nwestresults = nwest(g_data.y, g_data.X, floor( 1.3*T^(1/2) ) );   
        
        post_b(g_index,:)=repmat(nwestresults.beta,length(g_index),1);
    end
end
%postB = kron( post_b, ones(T,1) );
mu=zeros(N,1);
for i=1:N
yi=ds.y(ds.N == i );
Xi=ds.X(ds.N == i );
% yi(1,:)=[];
% Xi(1,:)=[];
mu(i,1)=mean(yi-Xi*post_b(i,:)');
%uNT_post((i-1)*T+1:i*T,:)=yi-Xi*post_b(i,:)'-mu(i,1);
end    
%post_Q =  mean(uNT_post.^2);
%% store
% result.group         = group;
result.est_post = post_b;
result.est_mu = mu;
result.est_classo = b;
% % result.est_lasso     = est_lasso;
% % % result.est_lasso_std = est_lasso_std;
% result.est_post_Q      = post_Q;
% % result.est_post_std  = est_post_std;
% % result.est_post_sig2 = est_post_sig2;
% result.est_lasso_Q   = classoQ;
% result.est_lasso_Q1   = classoQ1;
end
