function EST_Group_Coeff1(FILENAME, NEWFILENAME)

warning('off','all')
warning

%% load data
load(FILENAME);
load(NEWFILENAME);   % get result.T0hat
% T0hat    = result.T0hat;
T0hat    = result.T0tilde_Post;
thispreK    = 2;
thispostK    = 3;

Mhat     = size(T0hat,2);
if Mhat == 0    
    result.est = [];
    result.est_lambda = [];
else

%% construct data format
code=1:N;
code=kron(code,ones(T,1));
code=reshape(code,N*T,1);
period=[1:T];
period=period';
period=repmat(period,N);
period=period(:,1);

ds = dataset( code, period, y, X);
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};


%% data transformation: First-Diff
ds_o = dataset( code, period, y, X );
ds_o.Properties.VarNames = {'N'  'T'  'y'  'X'};

for i = 1:N
    yi = y(ds_o.N == i);
    Xi = X(ds_o.N == i, : );

    % first difference
    yid = [nan; yi(2:T,1) - yi(1:T-1,1)];
    Xid = [nan(1,p); Xi(2:T,:) - Xi(1:T-1,:)];

    y(ds_o.N==i) = yid;
    X(ds_o.N == i, :) =Xid;
end
%
% remove the 1st nan element
y(isnan(y)) = [];
X(isnan(X)) = [];


% % reset code and period
code=1:N;
code=kron(code,ones(T-1,1));
code=reshape(code,N*(T-1),1);
period=[2:T];
period=period';
period=repmat(period,N);
period=period(:,1);
%
% % new dataset
ds_o = dataset( code, period, y, X );
ds_o.Properties.VarNames = {'N'  'T'  'y'  'X'};

%% estimation given split points T0est
T0hat000 = [1 T0hat T];     % for original data
dT=3;
for thism = 1 : Mhat+1  % adjust for "thism"-th break 
    % thissample: T0hat^{thism-1} : T0hat^{thism+1}
    s1 = T0hat000(thism)+dT;
    s2 = T0hat000(thism+1)-dT;
    thisT = s2-s1+1;

    % subsample s1:s2
    thisds = ds( ds.T >= (s1) & ds.T<s2, : );
    thisds_o = ds_o( ds_o.T >= (s1) & ds_o.T<s2, : );
    % thisresult = est_Group2(thisds, Kmax);        
    
    % new lambda
    lambda0    = thisT^(-1/3);  % old lambda
    if s2<T0hat
    thisresult = est_Group2_lambda(thisds_o, thisds, thispreK, lambda0);
    else
    thisresult = est_Group2_lambda(thisds_o, thisds, thispostK, lambda0);
    end
    % save
    result.est(thism)           = thisresult;
end


end  




%% Save result
save(NEWFILENAME,'result')    

end

