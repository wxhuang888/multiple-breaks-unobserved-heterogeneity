%% 1. Setup and Environment
warning('off', 'all');  % Suppress all warnings
clear;                  % Clear workspace variables

%% 2. Load Data
FILENAME = 'sample10_all.mat'; % Specify the data file
load(FILENAME);

% Extract relevant variables
monthdate = data(:, 3);
min_val = min(monthdate);
max_val = max(monthdate);
Tmax = max_val - min_val + 1;
X = data(:, 5:8);  % Independent variables
y = data(:, 4)/100;     % Dependent variable, robust to y/100 scale
ID = data(:, 2);    % Individual IDs

%% 3. Summary statistics
datause = [data(:, 4) data(:, 5) data(:, 6) data(:, 7) data(:, 8)*100];

[n, K] = size(datause);

% Initialize arrays to store statistics
means = zeros(1, K);
stds = zeros(1, K);
p5 = zeros(1, K);
p25 = zeros(1, K);
medians = zeros(1, K);
p75 = zeros(1, K);
p95 = zeros(1, K);
%counts = zeros(1, K);

% Calculate statistics for each variable
for k = 1:K
    % Count non-NaN observations
    %counts(k) = sum(~isnan(datause(:,k)));

    means(k) = mean(datause(:,k), 'omitnan');
    stds(k) = std(datause(:,k), 'omitnan');
    p5(k) = prctile(datause(:,k), 5);
    p25(k) = prctile(datause(:,k), 25);
    medians(k) = median(datause(:,k), 'omitnan');
    p75(k) = prctile(datause(:,k), 75);
    p95(k) = prctile(datause(:,k), 95);
end

% Create variable names (V1, V2, ...)
varNames = cell(1, K);
for k = 1:K
    varNames{k} = ['V', num2str(k)];
end

% Create table with summary statistics
summary_stats = table(means', stds', p5', p25', medians', p75', p95', ...
    'RowNames', varNames, ...
    'VariableNames', {'Mean', 'SD', 'P5', 'P25', 'Median', 'P75', 'P95'});

% Display the summary statistics
disp('Summary Statistics:');
disp(summary_stats);


