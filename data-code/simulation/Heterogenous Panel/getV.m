function V = getV(ds, thiss1, thiss2, thisT2, thisdT, thislambda0)

thisT = thiss2 - thiss1 + 1;
thisN = size( unique(ds.N), 1);
%thisp = size( ds.X,2 );  % dimension of parameters

V = nan(thiss2,1);  %V0 = nan(thiss2,1); V2 = nan(thiss2,1);
for thistt = thiss1 + thisT2 : thisdT : thiss2 - thisT2   % searching range
    % thistt
    beta_pre = [];
    beta_post = [];
    res_pre = [];
    res_post = [];
    for i = 1:thisN
        % Data
        g_data_pre_i = ds( ds.N == i & ds.T >= thiss1 & ds.T < thistt, : );
        g_data_post_i = ds( ds.N == i & ds.T >= thistt & ds.T <= thiss2, : );

        % Estimation of \tilde{\beta}^{a/b}
        beta_pre(:,i)  = regress( g_data_pre_i.y, g_data_pre_i.X );
        % beta_pre(:,i)  = lsqlin(g_data_pre_i.X, g_data_pre_i.y, [], [], [], [], -log(thisT) * ones(size(g_data_pre_i.X,2),1), log(thisT)* ones(size(g_data_pre_i.X,2),1) );
        res_pre(:,i)   = g_data_pre_i.y - g_data_pre_i.X * beta_pre(:,i);

        beta_post(:,i) = regress( g_data_post_i.y, g_data_post_i.X );
        %beta_post(:,i) = lsqlin(g_data_post_i.X, g_data_post_i.y, [], [], [], [], -log(thisT) * ones(size(g_data_post_i.X,2),1), log(thisT)* ones(size(g_data_post_i.X,2),1) );
        res_post(:,i)   = g_data_post_i.y - g_data_post_i.X * beta_post(:,i);
    end

    Q_pre     = nansum(nansum(res_pre.^2));
    Q_post    = nansum(nansum(res_post.^2));

    % Calculate V
    V(thistt,1) = Q_pre + Q_post;
    if thistt ~= thiss1 + thisT2  & thistt ~= thiss2 - thisT2 
        V(thistt,1) = V(thistt,1) + thislambda0;
    end

end


end
