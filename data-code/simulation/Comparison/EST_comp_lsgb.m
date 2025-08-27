function EST_comp_lsgb(FILENAME, NEWFILENAME)
% This function uses Lumsdaine, Okui, and Wang (2023) estimation procedure

warning('off','all')
warning

%% load data
load(FILENAME);

% use FD transformation
K                       = size(X,2);
for k = 1:K
    Xmat(:,:,k)            = reshape(X(:,k), T, N)';
end
ymat                       = reshape(y,T,N)';

% FD transformation 
ymat = ymat(:,2:end) - ymat(:,1:end-1);
Xmat = Xmat(:,2:end,:) - Xmat(:,1:end-1,:); 

%% sequential
GB                      = 2;                                                % number of groups before the break
GA                      = 2;                                                % number of groups after the break
B                       = 1;                                                % number of breaks
[k_lsbg,gpB_lsbg,gpA_lsbg,bB_lsbg,bA_lsbg,seB_lsbg,seA_lsbg,~]   = est_gpbk_fd_tvG(ymat,Xmat,B,GB,GA);
% k_lsbg = estimated break point
% gpB_lsbg = estimated group memberships before the break
% gpA_lsbg = estimated group memberships after the break
% bB_lsbg = estimated slope coefficients before the break
% bA_lsbg = estimated slope coefficients after the break
% seB_lsbg = estimated standard deviation before the break
% seA_lsbg = estimated standard deviation after the break

result.T0hat    = k_lsbg;

%% Save result 

save(NEWFILENAME,'result')    
    

end

