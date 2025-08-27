function EST_Hetero_Inst_rolling(FILENAME, NEWFILENAME, est)

warning('off','all')
warning

%% load data
load(FILENAME);

%% tunings
r1      = floor( est.c * T^( est.c1 ) ) ;
lambda0 = est.clambda0 * N^( est.cN ) * T^( est.cT ) ; 
dT      = 1;   

%% construct data format
code=1:N;
code=kron(code,ones(T,1));
code=reshape(code,N*T,1);
period=[1:T];
period=period';
period=repmat(period,N);
period=period(:,1);

CX = [X ones(N*T,1)];  

ds = dataset( code, period, y, CX);
ds.Properties.VarNames = {'N'  'T'  'y'  'X'};

%% SEARCH T0hat
% The kth window is [(k − 1)r1, (k + 1)r1]
% k∈{1,...,⌊T/r1⌋−1}
% within each window [thiss1, thiss2]
% search thistt = thiss1 + 1 : Delta T : thiss2
dd = 0;   % check strictly extremes

% Start timing
tic;

% To store the estimate break dates
T0hat_long = [];
T0hat_short = [];

% Number of windows
num_windows = floor(T / r1);

% Whether breaks detected in each window
capt_inst = zeros(1, num_windows);

% To store the windows corresponding to T0hat_short
windowU = []; % Initialize an empty array to store windows

% Loop through each window
for k = 1:num_windows

    % Define the window range [thiss1, thiss2]
    thiss1 = (k - 1) * r1 + 1; 
    thiss2 = min((k + 1) * r1, T); 
    window = [thiss1, thiss2];

    % Tuning parameters
    thisT2 = size(ds.X, 2); % Number of regressors
    thisdT = dT;            % Step size
    thislambda0 = lambda0;  % Regularization parameter
    
    % Get V (objective function values for detecting breaks)
    thisV = getV(ds, thiss1, thiss2, thisT2, thisdT, thislambda0);

    % Check whether there is a break (T0hat)
    thisT0hat = find(thisV == min(thisV)); % Find the index of the minimum value in V
    if isempty(thisT0hat)
        thisT0hat = nan; % No break detected
        continue;        % Skip the rest of the loop for this window
    end

    % Check if the detected break is within the valid range
    if min(abs(thisT0hat - thiss1 - thisT2), abs(thiss2 - thisT2 - thisT0hat)) > dd
        % Store results if a valid break is detected
        capt_inst(k) = 1; % Mark this window as having a break
        T0hat_long = [T0hat_long, thisT0hat];
        T0hat_short = [T0hat_short, thisT0hat];
        windowU = [windowU; window]; % Append the window to windowU
    end 

    % Drop Method 1: Check for overlapping breaks
    if k > 1 && capt_inst(k) == 1 && capt_inst(k-1) == 1 && ...
       T0hat_long(end) < thiss1 + r1 && T0hat_long(end-1) > thiss1 - 2*r1
        % Remove the last detected break if it is redundant
        T0hat_short = T0hat_short(1:end-1);
        windowU(end, :) = []; % Remove the corresponding window from windowU
    end
end

% End timing
elapsedTime = toc; % Store the elapsed time


%% store the results
result.T0hat_long = T0hat_long;     
result.T0hat = T0hat_short;    
result.T0 = regime0;
result.WindowU = windowU;
result.r1 = r1;   
result.Time_Detect = elapsedTime; 

%% Save result 
save(NEWFILENAME,'result')    
    

end

