function EST_comp_gagfl(FILENAME, NEWFILENAME, G)
% This function uses Okui and Wang 2021 estimation procedure
% addpath('reference/OkuiWang2021replicationcode/replication files/simulation_replication')

warning('off','all')
warning

%% load data
load(FILENAME);

code=1:N;
code=kron(code,ones(T,1));
code=reshape(code,N*T,1);
period=[1:T];
period=period';
period=repmat(period,N);
period=period(:,1);

CX = X; % No constant in group
% CX = [X ones(N*T,1)];  % add constant
ds = dataset( code, period, y, CX);
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

%% sequential
option.maxLambda            = 100;
option.minLambda            = 0.01;
option.nGrid                = 40;

% %--- pls estimation (common breaks)
% b_init                                 = post_pls(y,CX,N,(1:T+1)');      % initial estimate of beta
% [regime_pls,alpha_pls,se_pls,~,~,~,~]  = pls(y,CX,N,b_init,option);

%--- gpls estimation
est_regime                                                      = [];
while isempty(est_regime)
    [beta_init,gi_init]                                         = gfe_est(ds,G);
    [est_regime,est_alpha,est_se,est_group,resQ,~]              = GPLS_est(ds,G,beta_init,gi_init,option); % grouped penalized least square (heterogeneous breaks)
end


% est_regime(1)   = [];
% est_regime(end) = [];
result.T0hat      = est_regime;
result.est_alpha  = est_alpha;
result.est_se     = est_se;
result.est_group  = est_group;

%% Save result 

save(NEWFILENAME,'result')    
    

end

