function result = est_Oracle(thisds, thisK, group0)
% this function get Q for c-lasso estimation

K = thisK;
N = max(thisds.N) - min(thisds.N) + 1;
T = max(thisds.T) - min(thisds.T) + 1;      
        
% prepare some parameters
y = thisds.y;
X = thisds.X;
p = size(X,2);

%% oracle estimation
NN = 1:N;

est_oracle     = zeros(K,p);
est_oracle_std = zeros(K,p);
group = logical(group0);

% 
for i = 1:K
    this_group= group(:,i);
    g_index = NN(this_group);

    if isempty(g_index)
        est_oracle(g_index,:)   =  nan;
        est_oracle_std(g_index,:)   =  nan;
    else
        g_data = thisds( ismember(thisds.N, g_index), : ); % group-specific data

        %oracle
        nwestresults = nwest(g_data.y, g_data.X, floor( 1.3*T^(1/2) ) );   
        est_oracle(i,:)=nwestresults.beta;
        est_oracle_std(i,:)=nwestresults.se;

    end

end

%% store

result.est_oracle     = est_oracle;
result.est_oracle_std = est_oracle_std;


end
