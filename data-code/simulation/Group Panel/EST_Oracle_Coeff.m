function EST_Oracle_Coeff(FILENAME, NEWFILENAME,Type)

warning('off','all')
warning

%% load data
load(FILENAME);
load(NEWFILENAME);   
% T0hat    = result.T0hat;

%% Define the true group membership and M/break dates
Mhat=1;
T0 = regime0(2);
% % %
if Type == 1
    thisKpre    = 2;
    thisKpost   = 2;
    N1A = 0.3*N;
    N1B = N1A;
    groupA=zeros(N,2);
    groupA(1:N1A,1)=1;
    groupA(N1A+1:N,2)=1;
    groupB=zeros(N,2);
    groupB(1:N1B,1)=1;
    groupB(N1B+1:N,2)=1;
elseif Type ==2
    thisKpre    = 2;
    thisKpost   = 2;
    N1A = 0.3*N;
    N1B = 0.5*N;
    groupA=zeros(N,2);
    groupA(1:N1A,1)=1;
    groupA(N1A+1:N,2)=1;
    groupB=zeros(N,2);
    groupB(1:N1B,1)=1;
    groupB(N1B+1:N,2)=1;
elseif Type ==3
    thisKpre    = 2;
    thisKpost   = 2;
    N1A = 0.3*N;
    N1B = 0.7*N;
    groupA=zeros(N,2);
    groupA(1:N1A,1)=1;
    groupA(N1A+1:N,2)=1;
    groupB=zeros(N,2);
    groupB(1:N1B,1)=1;
    groupB(N1B+1:N,2)=1;
elseif Type ==4
    thisKpre    = 2;
    thisKpost   = 3;
    N1A = 0.4*N;
    N1B = 0.4*N;
    N2B = 0.7*N;
    groupA=zeros(N,2);
    groupA(1:N1A,1)=1;
    groupA(N1A+1:N,2)=1;
    groupB=zeros(N,3);
    groupB(1:N1B,1)=1;
    groupB(N1B+1:N2B,2)=1;
    groupB(N2B+1:N,3)=1;
end
%% construct data format
code=1:N;
code=kron(code,ones(T,1));
code=reshape(code,N*T,1);
period=[1:T];
period=period';
period=repmat(period,N);
period=period(:,1);

ds_o = dataset(code, period, y, X);
ds_o.Properties.VarNames = {'N' 'T' 'y' 'X'};
for i = 1:N
    yi = y(ds_o.N == i);
    Xi = X(ds_o.N == i, :);
    
    % calculate means
    yi_mean = mean(yi);
    Xi_mean = mean(Xi);
    
    % demean transformation
    yi_demean = yi - yi_mean;
    Xi_demean = Xi - Xi_mean;
    
    % update the data
    y(ds_o.N == i) = yi_demean;
    X(ds_o.N == i, :) = Xi_demean;
end

ds_o = dataset(code, period, y, X);
ds_o.Properties.VarNames = {'N' 'T' 'y' 'X'};

%% estimation given split points T0est
T0hat000 = [1 T0 T];     % for original data
dT=3;
for thism = 1 : Mhat+1  % adjust for "thism"-th break
    % thissample: T0hat^{thism-1} : T0hat^{thism+1}
    s1 = T0hat000(thism)+dT;
    s2 = T0hat000(thism+1)-dT;
    %thisT = s2-s1+1;

    % subsample s1:s2
    thisds = ds_o( ds_o.T >= s1 & ds_o.T<s2, : );
    if s2<T0
        group0=groupA;
        thisresult = est_Oracle(thisds, thisKpre, group0);
    else
        group0=groupB;
        thisresult = est_Oracle(thisds, thisKpost, group0);
    end

    % save
    result.est_oracle(thism)           = thisresult;
end



%% Save result
save(NEWFILENAME,'result')

end

