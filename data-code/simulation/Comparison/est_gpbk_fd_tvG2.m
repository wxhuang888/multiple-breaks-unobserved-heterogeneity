function[k,gpmemB,gpmemA,coefB,coefA,seB,seA,resQB, resQA] = est_gpbk_fd_tvG2(ymat,Xmat,B,GB,GA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method simultaneously estimate the group structure, break point, and
% slope coefficients for panel models with fixed effects. 
% It allows that the number of groups varies over time.
%
% Input: 
% ymat: dependent variable in matrix form
% Xmat: explanatory variables in matrix form
% B: number of breaking points
%
% Output:
% k: estimated break point
% gpmemB: group membership estimate before the break
% gpmemA: group membership estimate after the break
% coefB: coefficient estimates before the break
% coefA: coefficient estimates after the break
% seB: standard error before the break
% seA: standard error after the break
% resQ: sum of squared residuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameterization
T       = size(ymat,2);
N       = size(ymat,1);
p       = size(Xmat,3);
S       = 30;                                                               % Maximum number of iterations for group estimation
SSR     = zeros(1,(T-2));
SSR(1)  = inf;
ymatfd  = diff(ymat,1,2);
Xmatfd  = diff(Xmat,1,2);

%% Break point estimation
for j = 1:B
    for bk = 2:(T-2)
        regime = [0,bk,T-1];                                                % Get the bounds of potential subperiods      
        for aa = 1:(j+1)                             
            T_frac          = length((regime(aa)+1):(regime(aa+1)));            
            y_frac          = reshape(ymatfd(:,(regime(aa)+1):(regime(aa+1)))',N*T_frac,1);            
            for l = 1:p            
                X_frac(:,l) = reshape(Xmatfd(:,(regime(aa)+1):(regime(aa+1)),l)',N*T_frac,1);                
            end            
            time_frac       = repmat(1:T_frac,1,N)';            
            code_frac       = kron((1:N)',ones(T_frac,1));                  
            idxdata_frac    = dataset(code_frac, time_frac, y_frac, X_frac);                        
            idxdata_frac.Properties.VarNames = {'N'  'T'  'y'  'X'};
            if aa == 1
                [~,~,ssr,]  = est_group(idxdata_frac,GB);
            elseif aa == 2
                [~,~,ssr,]  = est_group(idxdata_frac,GA);
            end                
            SSR(bk)         = SSR(bk) + ssr;        
            clear X_frac y_frac idxdata_frac time_frac code_frac
        end
    end
    [resQ,k]       = min(SSR);                                              % Compare SSR to get a breakdate
end
regime          = [0,k,T-1];
coefB           = zeros(GB,p);
coefA           = zeros(GA,p);
seB             = zeros(GB,p);
seA             = zeros(GA,p);
gpmem           = [];

%% Post break estimation: Classification for each subperiod
for j = 1:(B+1)            
    for aa = 1:j                                 
        T_frac          = length((regime(aa)+1):(regime(aa+1)));                    
        y_frac          = reshape(ymatfd(:,(regime(aa)+1):(regime(aa+1)))',N*T_frac,1);                   
        for l = 1:p                    
            X_frac(:,l) = reshape(Xmatfd(:,(regime(aa)+1):(regime(aa+1)),l)',N*T_frac,1);                            
        end        
        time_frac       = repmat(1:T_frac,1,N)';        
        code_frac       = kron((1:N)',ones(T_frac,1));                          
        idxdata_frac    = dataset(code_frac, time_frac, y_frac, X_frac);                                
        idxdata_frac.Properties.VarNames = {'N'  'T'  'y'  'X'};
        if aa == 1
            [coefB,gpmemB,resQB,seB]     = est_group(idxdata_frac,GB);
        elseif aa == 2
            [coefA,gpmemA,resQA,seA]     = est_group(idxdata_frac,GA);
        end
        clear X_frac y_frac idxdata_frac time_frac code_frac
    end
end