function EST_comp_gfl(FILENAME, NEWFILENAME)


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

% new dataset
ds = dataset( code, period, y, X);
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

%%
regime_fgl = nan(N,20);
for i = 1:N
    yi = y(ds.N == i);
    Xi = X(ds.N == i, : );
    Zi = ones(size(yi,1),1);

    %
    [regime,alpha,Sigma,R2,ssr,resid]=gfl(yi,Xi,Zi);
    % 
    regime_fgl(i,1:size(regime,1)-2) = regime(2:end-1,1)';
end

%% 
result.T0hat    = regime_fgl;
%% Save result 

save(NEWFILENAME,'result')    
    

end

