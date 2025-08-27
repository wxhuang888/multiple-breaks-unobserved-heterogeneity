function EST_Group_Refine_FERobust_NEW(FILENAME, NEWFILENAME)

warning('off','all')
warning

%% load data
load(FILENAME);
load(NEWFILENAME);   % get result.T0hat
T0hat    = result.T0hat;
Kmax = 2;

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

%% refine step
% Start timing
tic;

%
Mhat     = size(T0hat,2);
if Mhat == 0
    T0tilde_c = T0hat;
    T0tilde_p=T0hat;

end

windowU = result.WindowU;
WholeWindow = [1, T];

% Initialize the output for uncovered windows
uncovered = [];

% Start from the beginning of WholeWindow
current_start = WholeWindow(1);

% Check if windowU starts at the beginning of WholeWindow
if windowU(1, 1) == WholeWindow(1)
    % Add NaN to mark this special case
    uncovered = [uncovered; NaN, NaN];
end

% Loop through each window in windowU
for i = 1:size(windowU, 1)
    % Get the start and end of the current windowU range
    start_u = windowU(i, 1);
    end_u = windowU(i, 2);

    % If there's a gap before this windowU range, add it to uncovered
    if current_start < start_u
        uncovered = [uncovered; current_start, start_u - 1];
    end

    % Check for the special case where the gap (K) is less than r1
    if i > 1
        prev_end = windowU(i-1, 2); % Get the end of the previous interval
        K = start_u - prev_end;    % Calculate the gap
        r1 = result.r1;
        if K < r1
            % Add NaN to mark the special case
            uncovered = [uncovered; NaN, NaN];
        end
    end

    % Update the current start to the end of the current windowU range + 1
    current_start = max(current_start, end_u + 1);
end

% If there's a gap after the last windowU range, add it to uncovered
if current_start <= WholeWindow(2)
    uncovered = [uncovered; current_start, WholeWindow(2)];
end

% Check if windowU ends at the end of WholeWindow
if windowU(end, 2) == WholeWindow(2)
    % Add NaN to mark this special case
    uncovered = [uncovered; NaN, NaN];
end

% subsample C-Lasso estimation
thisbetaALL = [];
thisbetaALL0 = [];
thismuALL = [];

for thism = 1 : Mhat+1

    % Get the sub-window
    subWin = uncovered(thism, :);

    % Check if subWin contains NaN
    if ~any(isnan(subWin))
        % If subWin is valid, process the data
        thisds_o = ds_o(ds_o.T > subWin(1) - 1 & ds_o.T < subWin(2) + 1, :);  % first-diff
        thisds = ds(ds.T > subWin(1) - 1 & ds.T < subWin(2) + 1, :);

        % Call your estimation function
        thisresult = est_Group_new(thisds_o, thisds, Kmax);

        % Store the results
        thisbetaALL(:, thism) = thisresult.est_post;
        thisbetaALL0(:, thism) = thisresult.est_classo;
        thismuALL(:, thism) = thisresult.est_mu;
    else
        % If subWin is invalid, assign NaN for consistency
        thisbetaALL(:, thism) = NaN(N, 1);  % Assign NaN to the column
        thisbetaALL0(:, thism) = NaN(N, 1);
        thismuALL(:, thism) = NaN(N, 1);
    end
end
thismu = nanmean(thismuALL,2);

%T0hat000 = [1 T0hat T];
% refwin   = max(1, round(T^(0.3)) );
% refwin = max(1, round(T^(0.4)) );
for thism = 1 : Mhat  % adjust for "thism"-th break
    % refine
    % thissample: T0hat^{thism-1} : T0hat^{thism+1}
    % s1 = T0hat000(thism);
    % s2 = T0hat000(thism+2);

    QA = nan(T,1);
    QB = nan(T,1);
    QA_o = nan(T,1);
    QB_o = nan(T,1);
    Rewindow=result.WindowU(thism, :);
    subWin1 = uncovered(thism, :);
    subWin2 = uncovered(thism+1, :);
    if any(isnan(subWin1+subWin2))
        T0tilde_p(thism) = T0hat(thism);
        T0tilde_c(thism) = T0hat(thism);
    else
        for ss = Rewindow(1) : Rewindow(2)
            % pre-ss subsample
            thisdsA_ss = ds( ds.T >= Rewindow(1) & ds.T<ss+1, : );
            thisdsoA_ss = ds_o( ds_o.T >= Rewindow(1) & ds_o.T<ss+1, : );
            if isempty(thisdsA_ss)
                QA(ss,1) = 0;
                QA_o(ss,1) = 0;
            else
                QA(ss,1)= sum( (thisdsA_ss.y - kron(thismu,ones(size(unique(thisdsA_ss.T),1),1)) - thisdsA_ss.X .* kron(thisbetaALL(:,thism),ones(size(unique(thisdsA_ss.T),1),1))).^2, 'all');
                QA_o(ss,1)=sum((thisdsoA_ss.y  - thisdsoA_ss.X .* kron(thisbetaALL0(:,thism),ones(size(unique(thisdsoA_ss.T),1),1))).^2, 'all');
            end

            % post-ss subsample
            thisdsB_ss = ds( ds.T > ss & ds.T<=Rewindow(2), : );
            thisdsoB_ss = ds_o( ds_o.T > ss & ds_o.T<=Rewindow(2), : );
            if isempty(thisdsB_ss)
                QB(ss,1) = 0;
                QB_o(ss,1) = 0;
            else
                QB(ss,1)= sum( (thisdsB_ss.y - kron(thismu,ones(size(unique(thisdsB_ss.T),1),1)) - thisdsB_ss.X .* kron(thisbetaALL(:,thism+1),ones(size(unique(thisdsB_ss.T),1),1))).^2, 'all');
                QB_o(ss,1)=sum((thisdsoB_ss.y  - thisdsoB_ss.X .* kron(thisbetaALL0(:,thism+1),ones(size(unique(thisdsoB_ss.T),1),1))).^2, 'all');
            end

        end

        Q = QA + QB;
        Q_o = QA_o + QB_o;
        T0tilde_p(thism) = find( Q == min(Q) )+1;
        T0tilde_c(thism) = find( Q_o == min(Q_o) )+1;
    end
end

% End timing
elapsedTime = toc; % Store the elapsed time


% store
result.T0tilde_Post = T0tilde_p;
result.T0tilde_Classo = T0tilde_c;
result.Time_Refine = elapsedTime;


%% Save result
save(NEWFILENAME,'result')

end

