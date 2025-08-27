%% 1. Setup and Environment
warning('off', 'all');  
clear;                  
%% 2. Load Data
FILENAME = 'sample10_all.mat'; % Specify the data file
load(FILENAME);

% Extract relevant variables
monthdate = data(:, 3);
min_val = min(monthdate);
max_val = max(monthdate);
Tmax = max_val - min_val + 1;
X = data(:, 5:8);  % Independent variables
y = data(:, 4)/100; 
ID = data(:, 2);    % Individual IDs

%% 3. Define the Demeaned Dataset
ds_o = dataset(ID, monthdate, y, X);
ds_o.Properties.VarNames = {'N', 'T', 'y', 'X'};

% Get unique individuals (IDs)
unique_IDs = unique(ds_o.N);

% Initialize arrays for demeaned variables
demeaned_y = zeros(size(ds_o.y));
demeaned_X = zeros(size(ds_o.X));

% Loop over each unique ID to demean variables
for i = 1:length(unique_IDs)
    current_ID = unique_IDs(i);
    idx = ds_o.N == current_ID; 

    % Demean response variable (y)
    demeaned_y(idx) = ds_o.y(idx) - mean(ds_o.y(idx));

    % Demean independent variables (X)
    demeaned_X(idx, :) = ds_o.X(idx, :) - mean(ds_o.X(idx, :));
end

ds_o.y = demeaned_y;
ds_o.X = demeaned_X;

%% 4. Compute Initial Values for Beta 
p = size(ds_o.X, 2); 
beta_hat0 = zeros(length(unique_IDs), p); 
for i = 1:length(unique_IDs)
    yi = ds_o.y(ds_o.N == unique_IDs(i));  
    Xi = ds_o.X(ds_o.N == unique_IDs(i), :); 
    row_indices = find(ds_o.N == unique_IDs(i)); 

    if rank(Xi' * Xi) == size(Xi, 2) 
        beta_hat0(i, :) = (Xi' * Xi) \ (Xi' * yi); 
        y_hat = Xi * beta_hat0(i, :)'; 
        uhat(row_indices) = yi - y_hat; 
    else
        beta_hat0(i, :) = NaN; 
        uhat(row_indices) = NaN; 
    end
end
mse = (nanmean(uhat .* uhat)); 
% Remove individuals with NaN estimates in beta_hat0
valid_ids = all(~isnan(beta_hat0), 2);
beta_hat0 = beta_hat0(valid_ids, :);
ds_o = ds_o(ismember(ds_o.N, unique_IDs(valid_ids)), :);
ds = ds_o; 

%% 5. Parameter Tuning
est.c = 1; 
est.c1 = 0.6;
r1 = floor(est.c * Tmax^(est.c1)); 
dT = 1;  
dd = 0;  

%% 6. Break Point Estimation
% Initialize storage for estimated breakpoints
T0hat_long = [];
T0hat_short = [];

% Calculate number of rolling windows
num_windows = floor(Tmax / r1) - 1;
capt_inst = zeros(1, num_windows); 

% To store the windows corresponding to T0hat_short
windowU = []; 

% Loop through rolling windows
for k = 1:num_windows
    %fprintf('Processing window %d/%d...\n', k, num_windows);

    % Define window range [thiss1, thiss2]
    thiss1 = (k - 1) * r1 + min_val;
    thiss2 = min((k + 1) * r1 + min_val, Tmax + min_val + 1);
    window = [thiss1, thiss2];

    g_data = ds(ds.T >= thiss1 & ds.T <= thiss2, :);
    thisN = size(unique(g_data.N), 1);

    thisT2 = 0;   
    thisdT = dT;  

    est.cN = 1;
    est.cT = 0.15;
    est.clambda0 = 1; 
    thislambda0 = est.clambda0 * thisN^(est.cN) * Tmax^(est.cT);

    % Compute break detection statistic
    thisV = getV(g_data, thiss1, thiss2, thisT2, thisdT, thislambda0);

    thisT0hat = find(thisV == min(thisV));
    if isempty(thisT0hat)
        thisT0hat = nan;
    end

    if min(abs(thisT0hat - thiss1 - thisT2), abs(thiss2 - thisT2 - thisT0hat)) > dd
        capt_inst(k) = 1;
        T0hat_long = [T0hat_long, thisT0hat];
        T0hat_short = [T0hat_short, thisT0hat];
        windowU = [windowU; window];
    end

    if k > 1 && capt_inst(k) == 1 && capt_inst(k-1) == 1 && ...
        T0hat_long(end) < thiss1 + r1 && T0hat_long(end-1) > thiss1 - 2 * r1
        T0hat_short = T0hat_short(1:end-1);
        windowU(end, :) = [];
    end

    % Drop second break if windows intersect
    if size(windowU, 1) > 1  
        prev_window = windowU(end-1, :); 
        curr_window = windowU(end, :);   
    
        if curr_window(1) < prev_window(2) && curr_window(2) > prev_window(1)
            T0hat_short = T0hat_short(1:end-1);
            windowU(end, :) = [];
        end
    end

end

%% 7. Store the Results
result.T0hat_long = T0hat_long;
result.T0hat = T0hat_short;
result.WindowU=windowU;
save('sample10_T0hat.mat', 'result')
%% 8. Group Refine breaks
load("sample10_T0hat.mat")
T0hat=result.T0hat; 
N=length(unique(ds.N));
Kmax=3;
Mhat = numel(T0hat);
% subsample C-Lasso estimation
thisbetaALL = NaN(N, 4, Mhat+1);  
thisbetaALL0 = NaN(N, 4, Mhat+1);
thismuALL = NaN(N, Mhat+1);
T0hat000 = [1, T0hat, max(ds.T)];

for thism = 1 : Mhat+1
    % Define a new valid sub-window based on segmenting `T0hat000`
    subWin2 = [T0hat000(thism), T0hat000(thism+1)];

    % Extract data for this new valid sub-window
    thisds_o = ds_o(ds_o.T >= subWin2(1) & ds_o.T <= subWin2(2), :);  
    thisds = ds(ds.T >= subWin2(1) & ds.T <= subWin2(2), :);

    % Call your estimation function
    thisresult = est_Group_new_emp(thisds_o, thisds, Kmax, mse);

    % Store the results
    thisbetaALL(:,:, thism) = thisresult.est_post;
    thisbetaALL0(:,:, thism) = thisresult.est_classo;
    thismuALL(:, thism) = thisresult.est_mu;

end
thismu = nanmean(thismuALL,2);

for thism = 1 : Mhat  % adjust for "thism"-th break

    QA = nan(Tmax,1);
    QB = nan(Tmax,1);
    QA_o = nan(Tmax,1);
    QB_o = nan(Tmax,1);
    Rewindow=result.WindowU(thism, :);

    for ss = Rewindow(1) : Rewindow(2)
        % pre-ss subsample
        thisdsA_ss = ds( ds.T >= Rewindow(1) & ds.T<ss+1, : );
        thisdsoA_ss = ds_o( ds_o.T >= Rewindow(1) & ds_o.T<ss+1, : );
        if isempty(thisdsA_ss)
            QA(ss-min_val+1,1) = 0;
            QA_o(ss-min_val+1,1) = 0;
        else
            QA(ss-min_val+1,1)= sum( (thisdsA_ss.y - kron(thismu,ones(size(unique(thisdsA_ss.T),1),1)) - thisdsA_ss.X .* kron(thisbetaALL(:,:,thism),ones(size(unique(thisdsA_ss.T),1),1))).^2, 'all');
            QA_o(ss-min_val+1,1)=sum((thisdsoA_ss.y  - thisdsoA_ss.X .* kron(thisbetaALL0(:,:,thism),ones(size(unique(thisdsoA_ss.T),1),1))).^2, 'all');
        end

        % post-ss subsample
        thisdsB_ss = ds( ds.T > ss & ds.T<=Rewindow(2), : );
        thisdsoB_ss = ds_o( ds_o.T > ss & ds_o.T<=Rewindow(2), : );
        if isempty(thisdsB_ss)
            QB(ss-min_val+1,1) = 0;
            QB_o(ss-min_val+1,1) = 0;
        else
            QB(ss-min_val+1,1)= sum( (thisdsB_ss.y - kron(thismu,ones(size(unique(thisdsB_ss.T),1),1)) - thisdsB_ss.X .* kron(thisbetaALL(:,:,thism+1),ones(size(unique(thisdsB_ss.T),1),1))).^2, 'all');
            QB_o(ss-min_val+1,1)=sum((thisdsoB_ss.y  - thisdsoB_ss.X .* kron(thisbetaALL0(:,:,thism+1),ones(size(unique(thisdsoB_ss.T),1),1))).^2, 'all');
        end

        Q = QA + QB;
        Q_o = QA_o + QB_o;

    end
    T0tilde_p(thism) = find( Q == min(Q) )+1+min_val;
    T0tilde_c(thism) = find( Q_o == min(Q_o) )+1+min_val;
end

% store
result.T0tilde_Post = T0tilde_p;
result.T0tilde_Classo = T0tilde_c;

save('sample10_T0hat.mat', 'result')

%% 9. Estimation Given Breakpoints
FILENAME1 = 'sample10_T0hat.mat';
load(FILENAME1);
T0hat = result.T0tilde_Post;

Kmax = 6;
T0hat000 = [min(ds_o.T), T0hat, max(ds_o.T)];
IC_total_post = NaN(Kmax, length(T0hat) + 1);
IC_total_classo = NaN(Kmax, length(T0hat) + 1);

% Loop over partitions
for thism = 1:length(T0hat) + 1
    s1 = T0hat000(thism);
    s2 = T0hat000(thism + 1);

    % Extract subsample
    thisds_o = ds_o(ds_o.T >= s1 & ds_o.T < s2, :);

    for k = 1:Kmax
        try
            thisresult = est_Group_new_emp(thisds_o, thisds_o, k,mse);
            IC_total_post(k, thism) = thisresult.est_post_Q;
            IC_total_classo(k, thism) = thisresult.est_lasso_Q;
        catch
            if k > 1
                IC_total_post(k, thism) = IC_total_post(k-1, thism);
                IC_total_classo(k, thism) = IC_total_classo(k-1, thism);
            end
            fprintf('Error in est_Group_new_emp at thism=%d, k=%d. Using previous results.\n', thism, k);
        end
    end
end

%% 10. Compute Information Criterion (IC)
N = length(unique(ds_o.N));
T = Tmax;
pen = 2/3 * (N*T)^(-1/2) * (1:Kmax)';
IC_final_post = log(IC_total_post) + pen;

%% 11. Final Estimation Given Optimal Split Points and Best K
Mhat = numel(T0hat);
regression_results = cell(1, Mhat + 1);

if Mhat > 0
    T0hat000 = [min(ds_o.T), T0hat, max(ds_o.T)];
     %for thism =Mhat+1
    for thism = 1:Mhat + 1
        s1 = T0hat000(thism);
        s2 = T0hat000(thism + 1);

        % Extract subsample
        thisds_o = ds_o(ds_o.T >= s1 & ds_o.T < s2, :);

        % Ensure IC_final_post has enough columns
        if size(IC_final_post, 2) >= thism
            [~, Kbest] = min(IC_final_post(:, thism));
        else
            Kbest = 1; % Default to 1 if IC_final_post index is out of range
        end

        % Estimate with Best K and Handle Errors
        try
            regression_results{thism} = est_Group_new_emp(thisds_o, thisds_o, Kbest, mse);
        catch ME
            fprintf('Error in est_Group_new_emp at segment %d. Skipping...\n', thism);
            regression_results{thism} = []; % Store empty if error occurs
        end
    end
end

save('sample10_est_groupIC.mat', 'IC_final_post')

save('sample10_est_group.mat', 'regression_results')

load('sample10_est_group.mat')
Re1_beta = unique(regression_results{1, 1}.est_post, 'rows', 'stable');
Re2_beta = unique(regression_results{1, 2}.est_post, 'rows', 'stable');
Re3_beta = unique(regression_results{1, 3}.est_post, 'rows', 'stable');
Re4_beta = unique(regression_results{1, 4}.est_post, 'rows', 'stable');

beta1 = [100*Re1_beta(1,:)' 100*Re1_beta(2,:)' 100*Re1_beta(3,:)'];
beta2 = [100*Re2_beta(1,:)' 100*Re2_beta(2,:)' 100*Re2_beta(3,:)' 100*Re2_beta(4,:)'];
beta3 = [100*Re3_beta(1,:)' 100*Re3_beta(2,:)' 100*Re3_beta(3,:)'];
beta4 = [100*Re4_beta(1,:)' 100*Re4_beta(2,:)' 100*Re4_beta(3,:)'];

Re1_beta_std = unique(regression_results{1, 1}.est_post_std1, 'rows', 'stable');
Re2_beta_std = unique(regression_results{1, 2}.est_post_std1, 'rows', 'stable');
Re3_beta_std = unique(regression_results{1, 3}.est_post_std1, 'rows', 'stable');
Re4_beta_std = unique(regression_results{1, 4}.est_post_std1, 'rows', 'stable');

beta1_std = [100*Re1_beta_std(1,:)' 100*Re1_beta_std(2,:)' 100*Re1_beta_std(3,:)'];
beta2_std = [100*Re2_beta_std(1,:)' 100*Re2_beta_std(2,:)' 100*Re2_beta_std(3,:)' 100*Re2_beta_std(4,:)'];
beta3_std = [100*Re3_beta_std(1,:)' 100*Re3_beta_std(2,:)' 100*Re3_beta_std(3,:)'];
beta4_std = [100*Re4_beta_std(1,:)' 100*Re4_beta_std(2,:)' 100*Re4_beta_std(3,:)'];

results.beta1 = beta1;
results.beta2 = beta2;
results.beta3 = beta3;
results.beta4 = beta4;

results.beta1_std = beta1_std;
results.beta2_std = beta2_std;
results.beta3_std = beta3_std;
results.beta4_std = beta4_std;
save('sample10_est_groupcoeff.mat', 'results');
%% Figure 1: ﻿The S&P 500 index with the estimated break dates
clear;
load('SP500.mat'); 
sp500_time = SP500.caldt;
sp500_value = SP500.spindx;

% Define breakpoints as datetime (obtain the exact dates from dates.xlsx)
% T0hat=result.T0tilde_Post; [674,706,754]
breakpoints = [datetime(2016, 3, 1), datetime(2018, 11, 1), datetime(2022, 11, 1)];

% Plot the S&P 500 index
figure;
plot(sp500_time, sp500_value, 'b-', 'LineWidth', 1.5); 
hold on;

for i = 1:length(breakpoints)

    xline(breakpoints(i), 'r--', 'LineWidth', 1.2); 
    

    text(breakpoints(i), max(sp500_value) * 0.95, datestr(breakpoints(i), 'yyyy-mmm'), ...
        'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center');
end

% Add labels and title
xlabel('Year','FontSize', 21);
ylabel('S&P 500 Index','FontSize', 21);
grid on;
xlim([datetime(2014,1,1), datetime(2023,12,31)]);
xtickformat('yyyy');
set(gca, 'FontSize', 16);
legend('S&P 500 Index', 'Break dates', 'Location', 'Best');

saveas(gcf, 'sp500.png')


%% ﻿Figure 2: The sub-panel coefficient estimates of risk premiums
clear;
FILENAME1 = 'sample10_T0hat.mat';
load(FILENAME1);

FILENAME2 = 'sample10_all.mat';
load(FILENAME2);
monthdate=data(:,3);
min_val = min(monthdate);
max_val = max(monthdate);
Tmax=max_val-min_val+1;
X=data(:,5:8);
y=data(:,4);
ID = data(:, 2);    % Individual IDs
%% Define the Demeaned Dataset
ds_o = dataset(ID, monthdate, y, X);
ds_o.Properties.VarNames = {'N', 'T', 'y', 'X'};

% Get unique individuals (IDs)
unique_IDs = unique(ds_o.N);

% Initialize arrays for demeaned variables
demeaned_y = zeros(size(ds_o.y));
demeaned_X = zeros(size(ds_o.X));

% Loop over each unique ID to demean variables
for i = 1:length(unique_IDs)
    current_ID = unique_IDs(i);
    idx = ds_o.N == current_ID; % Identify relevant rows

    % Demean response variable (y)
    demeaned_y(idx) = ds_o.y(idx) - mean(ds_o.y(idx));

    % Demean independent variables (X)
    demeaned_X(idx, :) = ds_o.X(idx, :) - mean(ds_o.X(idx, :));
end

% Store demeaned variables back in the dataset
ds_o.y = demeaned_y;
ds_o.X = demeaned_X;

%% Compute Initial Values for Beta
p = size(ds_o.X, 2); % Number of predictors
beta_hat0 = zeros(length(unique_IDs), p);

for i = 1:length(unique_IDs)
    yi = ds_o.y(ds_o.N == unique_IDs(i));
    Xi = ds_o.X(ds_o.N == unique_IDs(i), :);

    if rank(Xi' * Xi) == size(Xi, 2) 
        beta_hat0(i, :) = (Xi' * Xi) \ (Xi' * yi); 
    else
        beta_hat0(i, :) = NaN; 
    end
end

% Remove individuals with NaN estimates in beta_hat0
valid_ids = all(~isnan(beta_hat0), 2);
beta_hat0 = beta_hat0(valid_ids, :);
ds_o = ds_o(ismember(ds_o.N, unique_IDs(valid_ids)), :);
ds = ds_o; % Store cleaned dataset

T0hat=result.T0tilde_Post;
%% Perform Subsample Panel Regression Using Demeaned Data
% Use breakpoints to segment the dataset
breakpoints = sort([min_val, T0hat, max_val + 1]); 
num_segments = length(breakpoints) - 1; 

% Initialize storage for regression results
regression_results_all = cell(1, num_segments);

for seg = 1:num_segments
    % Define the segment range
    segment_start = breakpoints(seg);
    segment_end = breakpoints(seg + 1) - 1; 
    segment_data = ds(ds.T >= segment_start & ds.T <= segment_end, :);
    y_segment = segment_data.y;
    X_segment = segment_data.X;

    if rank(X_segment' * X_segment) == size(X_segment, 2)
        beta_hat = (X_segment' * X_segment) \ (X_segment' * y_segment); 
        residuals = y_segment - X_segment * beta_hat; 
        RSS = sum(residuals.^2); 
        n = size(X_segment, 1); 
        p = size(X_segment, 2); 
        sigma_squared = RSS / (n - p); 
        var_beta = sigma_squared * inv(X_segment' * X_segment); 
        std_errors = sqrt(diag(var_beta)); 
        R2 = 1 - sum(residuals.^2) / sum((y_segment - mean(y_segment)).^2); 
    else
        beta_hat = nan(size(X_segment, 2), 1); 
        residuals = nan(size(y_segment));
        std_errors = nan(size(X_segment, 2), 1); 
        R2 = nan;
    end

    % Store regression results
    regression_results_all{seg}.beta_hat = beta_hat';
    regression_results_all{seg}.std_errors = std_errors';
    regression_results_all{seg}.R2 = R2;
    regression_results_all{seg}.segment_range = [segment_start, segment_end];
    regression_results_all{seg}.num_obs = size(X_segment, 1);
end

save('sample10_est_all.mat', 'regression_results_all')

FILENAME1 = 'sample10_est_all.mat';
load(FILENAME1);
beta = [
    regression_results_all{1, 1}.beta_hat;
    regression_results_all{1, 2}.beta_hat;
    regression_results_all{1, 3}.beta_hat;
    regression_results_all{1, 4}.beta_hat;
    ];

beta(:,2)=-beta(:,2);

std = [
    regression_results_all{1, 1}.std_errors;
    regression_results_all{1, 2}.std_errors;
    regression_results_all{1, 3}.std_errors;
    regression_results_all{1, 4}.std_errors;
    ];

load('sample10_est_group.mat')
Re1_beta=unique(regression_results{1, 1}.est_post,'rows');
Re2_beta=unique(regression_results{1, 2}.est_post,'rows');
Re3_beta=unique(regression_results{1, 3}.est_post,'rows');
Re4_beta=unique(regression_results{1, 4}.est_post,'rows');

beta1 = [
     100*Re1_beta(1,:);
     100*Re2_beta(1,:);
     100*Re3_beta(1,:);
     100*Re4_beta(1,:);
    ];
beta1(:,2)=-beta1(:,2);
beta2 = [
     100*Re1_beta(2,:);
     100*Re2_beta(2,:);
     100*Re3_beta(2,:);
     100*Re4_beta(2,:);
    ];
beta2(:,2)=-beta2(:,2);
beta3 = [
    100*Re1_beta(3,:);
    100*Re2_beta(3,:);
    100*Re3_beta(3,:);
    100*Re4_beta(3,:);
    ];
beta3(:,2)=-beta3(:,2);

beta4 = [
    NaN NaN NaN NaN;
    100*Re2_beta(4,:);
    NaN NaN NaN NaN;
    NaN NaN NaN NaN;
    ];
beta4(:,2)=-beta4(:,2);

breakpoints = { '201401-201602', '201603-201810', '201811-202210', ...
'202211-202312' };
time_ranges = [
   2014.0000   2016.1667
   2016.1667   2018.8333
   2018.8333   2022.8333
   2022.8333   2024.0000];

% Extract group sizes
Group1 = sum(regression_results{1, 1}.group); % [203, 847, 293]
Group2 = sum(regression_results{1, 2}.group); % Add your values here
Group3 = sum(regression_results{1, 3}.group); % Add your values here
Group4 = sum(regression_results{1, 4}.group); % Add your values here

% Store group sizes in a matrix
group_sizes = {Group1, Group2, Group3, Group4};

% Labels for the plots
titles = { '\bf{Market Beta}', '\bf{Size}', '\bf{Book to Market Ratio}', '\bf{Momentum}' };
y_labels = { '\bf{$\lambda_{Mkt}$}', '\bf{$-\lambda_{Size}$}', '\bf{$\lambda_{BM}$}', '\bf{$\lambda_{Mom}$}' };

% Initialize the figure
figure;

% Marker styles and grey color for all markers
marker_styles = {'o', 'o', 'o', 'o'}; 
marker_color = [0.5, 0.5, 0.5]; 

% Base marker size and scaling factor
base_marker_size = 3;
scale_factor = 0.01; 

% Loop over each coefficient to generate the plots
for i = 1:4
    subplot(2, 2, i); 
    x = []; % Time points
    y = []; % Beta coefficients
    for j = 1:size(time_ranges, 1)
        x = [x, time_ranges(j, 1), time_ranges(j, 2)]; 
        y = [y, beta(j, i), beta(j, i)]; 
    end

    
    plot(x, y, 'k-', 'LineWidth', 1.5); 
    hold on;

    for j = 1:size(time_ranges, 1)
        % Upper and lower bounds
        std_up = beta(j, i) + 1.96 * std(j, i); % Upper bound
        std_down = beta(j, i) - 1.96 * std(j, i); % Lower bound

        plot([time_ranges(j, 1), time_ranges(j, 2)], [std_up, std_up], 'r--', 'LineWidth', 1.2); % Upper bound
        plot([time_ranges(j, 1), time_ranges(j, 2)], [std_down, std_down], 'r--', 'LineWidth', 1.2); % Lower bound
    end

    for k = 1:4
        beta_current = eval(['beta' num2str(k)]);
        
        % Plot markers for the current beta dataset
        for j = 1:size(time_ranges, 1)
            midpoint = (time_ranges(j, 1) + time_ranges(j, 2)) / 2;           
            if ~isnan(beta_current(j, i))
                
                current_group_size = 0;
                if k <= length(group_sizes{j})
                    current_group_size = group_sizes{j}(k);
                end
                
                marker_size = base_marker_size + scale_factor * current_group_size;
                                
                plot(midpoint, beta_current(j, i), marker_styles{k}, 'MarkerEdgeColor', marker_color, ...
                    'MarkerFaceColor', marker_color, 'MarkerSize', marker_size);
            end
        end
    end

    xlabel('Year');
    ylabel(y_labels{i}, 'Interpreter', 'latex', 'FontSize', 18, 'FontWeight', 'bold');
    title(titles{i},'FontSize', 16);

    grid on;

    % Set x-axis limits to end at 2023m12
    xlim([2014, 2024]);
    
    set(gca, 'FontSize', 15); 
end

saveas(gcf, 'betaestimatesGroup.png')


%% ﻿Figure 12: ﻿The information criterion results across the four regimes
clear;
load("sample10_est_groupIC.mat")


% Titles for each subplot
titles = { ...
    '201401-201602', ...
    '201603-201810', ...
    '201811-202210', ...
    '202211-202312' ...
};

figure;

% Loop through each column and create a subplot
for col = 1:size(IC_final_post, 2)
    subplot(2, 2, col);
    plot(IC_final_post(:, col), 'LineWidth', 1.5); 
    title(titles{col}); 
    xlabel('Number of Groups'); 
    ylabel('IC Values'); 
    grid on; 
end

saveas(gcf, 'numGroups.png')