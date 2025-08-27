function EST_comp_pls(FILENAME, NEWFILENAME)


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

% CX = [ones(N*T,1) X]; 

% %% first-difference data
ds_o = dataset( code, period, y, X );
ds_o.Properties.VarNames = {'N'  'T'  'y'  'X'};

for i = 1:N
    yi = y(ds_o.N == i);
    Xi = X(ds_o.N == i, : );
  %  ui = u(ds_o.N == i);

    % first difference
    yid = [nan; yi(2:T,1) - yi(1:T-1,1)];
    Xid = [nan(1,p); Xi(2:T,:) - Xi(1:T-1,:)];
   % uid = [nan; ui(2:T,1) - ui(1:T-1,1)];

    y(ds_o.N==i) = yid;
    X(ds_o.N == i, :) =Xid;
  %  u(ds_o.N==i)= uid;
end

% remove the 1st nan element
y(isnan(y)) = [];
X(isnan(X)) = [];
%u(isnan(u)) = [];

% reset code and period
code=1:N;
code=kron(code,ones(T-1,1));
code=reshape(code,N*(T-1),1);
period=[2:T];
period=period';
period=repmat(period,N);
period=period(:,1);

% new dataset
ds = dataset( code, period, y, X);
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};
% z=[];
%% sequential
option.maxLambda            = 100;
option.minLambda            = 0.01;
option.nGrid                = 40;



%--- pls estimation (common breaks)
b_init                                 = post_pls(y,X,N,(1:T)');      % initial estimate of beta
[regime_pls,~,~,~,~,~,~]  = pls(y,X,N,b_init,option);

regime_pls(1)   = [];
regime_pls(end) = [];
result.T0hat    = regime_pls;
% weight=ones(T-1,1);
% XTol=1e-4;
% maxIter=100;
% lambda = 25;
% [tht1,~]=hpbcdq(y,X,N,z,lambda,weight,XTol,maxIter);
% theta = reshape(tht1(:,1),N,T-1)';
% [regime,~]=findbreaks(theta,[],2);
% regime(1)   = [];
% regime(end) = [];
% result.T0hat    = regime;
%% Save result 

save(NEWFILENAME,'result')    
    

end

