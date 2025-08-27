%% 1. Setup and Environment
warning('off', 'all');  % Suppress all warnings
clear;                  % Clear workspace variables
addpath(genpath('reference'));
%% 2. Load Data
FILENAME = 'sample10_all.mat'; % Specify the data file
load(FILENAME);

% Extract relevant variables
monthdate = data(:, 3);
min_val = min(monthdate);
max_val = max(monthdate);
Tmax = max_val - min_val + 1;
X = data(:, 5:8);  % Independent variables
y = data(:, 4)/100;     % Dependent variable
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
    idx = ds_o.N == current_ID; % Identify relevant rows

    % Demean response variable (y)
    demeaned_y(idx) = ds_o.y(idx) - mean(ds_o.y(idx));

    % Demean independent variables (X)
    demeaned_X(idx, :) = ds_o.X(idx, :) - mean(ds_o.X(idx, :));
end

% Store demeaned variables back in the dataset
ds_o.y = demeaned_y;
ds_o.X = demeaned_X;

%% 4. Compute Initial Values for Beta and Residuals (uhat)
p = size(ds_o.X, 2); 
beta_hat0 = zeros(length(unique_IDs), p); 
uhat = zeros(length(ds_o.y), 1); 

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
valid_ids = all(~isnan(beta_hat0), 2);
beta_hat0 = beta_hat0(valid_ids, :);
ds_o = ds_o(ismember(ds_o.N, unique_IDs(valid_ids)), :);
ds = ds_o; % Store cleaned dataset
unique_IDs = unique(ds.N);
[~, new_ID] = ismember(ds.N, unique_IDs);
ds.N = new_ID;

variableNames = {'WindowSize', 'T0hat_B1', 'T0hat_B2', 'T0hat_B3'};
resultsTable = table('Size', [0, 4], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames', variableNames);

load("sample10_T0hat.mat")
T0hat_M=result.T0tilde_Post;
for window = 5:12
    % Define the time intervals for each subsample
    T1_start = T0hat_M(1) - window; T1_end = T0hat_M(1) + window; % Subsample 1
    T2_start = T0hat_M(2) - window; T2_end = T0hat_M(2) + window; % Subsample 2
    T3_start = T0hat_M(3) - window; T3_end = T0hat_M(3) + window; % Subsample 3

    % --- Subsample 1 ---
    ds1 = ds(ds.T >= T1_start & ds.T <= T1_end, :);
    unique_ID1 = unique(ds1.N);
    [~, new_ID1] = ismember(ds1.N, unique_ID1);
    ds1.N = new_ID1;
    N1 = length(unique_ID1);
    T1 = T1_end - T1_start + 1;
    K = size(ds1.X,2);
    
    X1 = ds1.X;
    y1 = ds1.y;
    xmat1 = reshape(X1, T1, N1, K);
    ymat1 = reshape(y1, T1, N1);
    zmat1 = xmat1;
    ifull = 1;
    M=1;
    [est_break_pts1] = Revise_bfk_multi_break_detect_M(N1,T1,ymat1,xmat1,zmat1,ifull, M);
    T0hat_B1 = T1_start + est_break_pts1; 
    
    % --- Subsample 2 ---
    ds2 = ds(ds.T >= T2_start & ds.T <= T2_end, :);
    unique_ID2 = unique(ds2.N);
    [~, new_ID2] = ismember(ds2.N, unique_ID2);
    ds2.N = new_ID2;
    N2 = length(unique_ID2);
    T2 = T2_end - T2_start + 1;
    K = size(ds2.X,2);
    
    X2 = ds2.X;
    y2 = ds2.y;
    xmat2 = reshape(X2, T2, N2, K);
    ymat2 = reshape(y2, T2, N2);
    zmat2 = xmat2;
    ifull = 1;
    M=1;
    [est_break_pts2] = Revise_bfk_multi_break_detect_M(N2,T2,ymat2,xmat2,zmat2,ifull, M);
    T0hat_B2 = T2_start + est_break_pts2; 

    % --- Subsample 3 ---
    ds3 = ds(ds.T >= T3_start & ds.T <= T3_end, :);
    unique_ID3 = unique(ds3.N);
    [~, new_ID3] = ismember(ds3.N, unique_ID3);
    ds3.N = new_ID3;
    N3 = length(unique_ID3);
    T3 = T3_end - T3_start + 1;
    K = size(ds3.X,2);
    
    X3 = ds3.X;
    y3 = ds3.y;
    xmat3 = reshape(X3, T3, N3, K);
    ymat3 = reshape(y3, T3, N3);
    zmat3 = xmat3;
    ifull = 1;
    M=1;
    [est_break_pts3] = Revise_bfk_multi_break_detect_M(N3,T3,ymat3,xmat3,zmat3,ifull, M);
    T0hat_B3 = T3_start + est_break_pts3; 

    newRow = {window, T0hat_B1, T0hat_B2, T0hat_B3};
    resultsTable = [resultsTable; newRow];

end

save('sample10_BFK.mat', 'resultsTable')