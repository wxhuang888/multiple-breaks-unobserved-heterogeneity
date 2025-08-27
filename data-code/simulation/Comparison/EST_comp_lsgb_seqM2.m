function EST_comp_lsgb_seqM2(FILENAME, NEWFILENAME)
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
B                       = 1;                                              % number of breaks

%% 1st break
[k_lsbg,gpB_lsbg,gpA_lsbg,bB_lsbg,bA_lsbg,seB_lsbg,seA_lsbg, resQB, resQA]   = est_gpbk_fd_tvG2(ymat,Xmat,B,GB,GA);
% k_lsbg = estimated break point
% gpB_lsbg = estimated group memberships before the break
% gpA_lsbg = estimated group memberships after the break
% bB_lsbg = estimated slope coefficients before the break
% bA_lsbg = estimated slope coefficients after the break
% seB_lsbg = estimated standard deviation before the break
% seA_lsbg = estimated standard deviation after the break

%% 2nd break
% in subsample 1:k_lsbg
thisymat = ymat( :, 1:k_lsbg );
thisXmat = Xmat( :, 1:k_lsbg );
[k_lsbg_1,gpB_lsbg,gpA_lsbg,bB_lsbg,bA_lsbg,seB_lsbg,seA_lsbg, thisresQB, thisresQA]   = est_gpbk_fd_tvG2(thisymat,thisXmat,B,GB,GA);
SSR1 = (thisresQB + thisresQA) + resQA;

% in subsample k_lsbg+1:T
thisymat = ymat( :, k_lsbg+1:end );
thisXmat = Xmat( :, k_lsbg+1:end );
[k_lsbg_2,gpB_lsbg,gpA_lsbg,bB_lsbg,bA_lsbg,seB_lsbg,seA_lsbg, thisresQB, thisresQA]   = est_gpbk_fd_tvG2(thisymat,thisXmat,B,GB,GA);
SSR2 = resQB + (thisresQB + thisresQA);


% compare SSR and pick the 2nd break
if SSR1 <= SSR2
    T0hat = [k_lsbg_1; k_lsbg];
elseif SSR1 > SSR2
    T0hat = [k_lsbg; k_lsbg+k_lsbg_2];
end



%% save
result.T0hat    = T0hat;

%% Save result 

save(NEWFILENAME,'result')    
    

end

