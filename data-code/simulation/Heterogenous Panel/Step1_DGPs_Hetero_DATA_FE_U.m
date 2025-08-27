
% This code generate Heterogenous Panel with M breaks in slope coefficient \beta
% for DGPI 1--4

warning off
clc
clear all

%% Replicate 500
for ii=1:1
    for DGPnum = [1 ]             % DGPI: 1,2,3,4
        %% setting
        for M    = [2]                % # breaks; # segments = M+1
            for N    = [50 ]         
                for T    = [100 ]   
                    for NF   = [0.5]   % fraction of individuals to remain unchanged across regimes

                        pathhtc = '';
                        disp(ii)
                        
                        tic;
                        p    = 1;                  % # of regressors
                        
                        % generate the beta parameter
                        sigma_a = 1; % standard deviation for alpha generation
                        jmin = 1;    % minimum step size for alpha increments
                        jmax = 2;    % maximum step size for alpha increments
                        RandomB=randn(1);

                        regime0 = [1;floor(linspace(T/(M+1),M*T/(M+1),M))'+1+floor(RandomB);T+1]; % true regime
                        N1F = NF*N;
                        alpha0 = simu_genalpha_N1F_RS(sigma_a,jmin,jmax,p,N,M,N1F);
                        beta0 = alpha2beta2(alpha0, regime0); % convert alpha0 to beta0

                        % signal - noise
                        sigma_e = 1; % standard deviation of noise
                        X = randn(N * T, p); % generate random design matrix X (Normal distribution)

                        mu = randn(N, 1); % generate random intercepts for each individual
                        mu0 = kron(mu, ones(T, 1)); % repeat mu for each period T
                        mu0 = reshape(mu0, N * T, 1); % reshape into a column vector
                        X = mu0+X;
                        % DGP: generate the dependent variable y
                        if DGPnum ==1
                            U = randn(N * T, 1);
                            y = mu0 + sum(X .* beta0, 2) + sigma_e * U; % final response variable
                        elseif DGPnum ==2
                            U = randn(N * T, 1);
                            e = U;
                            for t = 2:T
                                for i = 1:N
                                    e((i-1)*T+t) = U((i-1)*T+t) + 0.2 * U((i-1)*T+t-1);
                                end
                            end
                            y = mu0+sum(X.*beta0,2) + sigma_e * e;
                        elseif DGPnum == 3
                            U = randn(N * T, 1);
                            e = U;
                            h = ones(N*T,1);
                            for t = 2:T
                                for i = 1:N
                                    h((i-1)*T+t) = 0.05 + 0.05 * (e((i-1)*T+t-1))^2 + 0.9 * h((i-1)*T+t-1);
                                    e((i-1)*T+t) = sqrt(h((i-1)*T+t)) * U((i-1)*T+t);
                                end
                            end
                            y = mu0+sum(X.*beta0,2) + sigma_e * e;
                        elseif DGPnum ==4
                            td = 8; % degrees of freedom for t-distribution
                            U = trnd(td, N * T, 1); % generate t-distribution noise
                            e = U;
                            for t = 2:T
                                for i = 1:N
                                    e((i-1)*T+t) = U((i-1)*T+t) + 0.2 * U((i-1)*T+t-1);
                                end
                            end
                            y = mu0+sum(X.*beta0,2) + sigma_e * e;
                        end
                        %% Save result
                        clear alpha0 beta0 mu0 mu U
                        DataFolder = fullfile( sprintf(append(pathhtc, 'DATA/DATA', num2str(DGPnum), '_U_M',num2str( M ),'_NF',num2str( NF ),'_N',num2str( N ),'_T', num2str( T ))) );
                        if ~exist(DataFolder, 'dir')
                            mkdir(DataFolder)
                        end

                        name = [DataFolder,'/DATA_',num2str(ii)];
                        save(name)

                    end
                end
            end

        end
    end
    save ii

    clear all

    load ii


    toc;
end







