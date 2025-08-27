
% This code generates Group Panel with M breaks in slope coefficient \beta
% for DGPG 1--4 
% for different arrival dates of the break T_1^0/T


warning off
clc
clear all

%% Replicate 500
for Type =[1 ]                 % DGPG 1--4
    for Bfrac=[0.1:0.1:0.9]    % T_1^0/T
        for N = [50 100]
            for T = [50 100]
                for ii=1:1
                    disp(ii)
                    
                    tic;
                    %% setting
                    RandomB=randn(1);
                    T1 = floor(Bfrac * T);      % break point
                    regime0 = [1;T1+1;T+1];
                    K = 2;    % # of groups
                    p = 1;    % # of regressors

                    if Type == 1
                        % only group coefficient change, no group member change
                        N1A = 0.3*N;
                        N1B = N1A;
                        N2 = N;   % only 2 groups
                        a0A = [-0.4;
                            1.6];
                        a0B = [-0.4;
                            0.4];



                    elseif Type == 2
                        % only group member change
                        N1A = 0.3*N;
                        N1B = 0.5*N;
                        N2 = N;   % only 2 groups

                        a0A = [-0.4;
                            1.6];
                        a0B = a0A;

                    elseif Type == 3
                        % mixture of Type 1 and Type 2
                        N1A = 0.3*N;
                        N1B = 0.7*N;
                        N2 = N;   % only 2 groups

                        a0A = [-0.4;
                            1.6];
                        % a0A = [-0.4;
                        %     0.4];
                        a0B = [-0.4;
                            0.8];

                    elseif Type == 4
                        % change in the number of groups
                        N1A = 0.4*N;
                        N1B = 0.4*N;
                        N2B = 0.7*N;
                        N2 = N;   %  3 group


                        a0A = [0.8;
                            1.6];
                        a0B = [-0.4;1.6;0.8];

                    elseif Type == 0
                        % no changes
                        N1A = 0.3*N;
                        N1B = N1A;
                        N2 = N;   % only 2 groups

                        a0A = [-0.4;
                            1.6];
                        a0B = a0A;

                    else
                        error('check Type of Instability')
                    end

                    % same as a0A, a0B, but for each series
                    beta0A = [ones(N1A,1) * a0A(1,:);
                        ones((N-N1A),1) * a0A(2,:)];

                    if Type ~= 4
                        beta0B = [ones(N1B,1) * a0B(1,:);
                            ones((N-N1B),1) * a0B(2,:)];
                    end
                    if Type == 4
                        beta0B = [ones(N1B,1) * a0B(1,:);
                            ones((N2B-N1B),1) * a0B(2,:);
                            ones((N2-N2B),1) * a0B(3,:)];
                    end

                    sigma_x = 1;
                    sigma_e = 1;

                    %% DGP
                    % Model: y_{it} =  x_{it}' \alpha_{kt} + u_{it}
                    %        \alpha_{kt} = \alpha_{k}^{A} if t < tau
                    %        \alpha_{kt} = \alpha_{k}^{B} if t < tau

                    % intercept
                    mu0=randn(N,1);
                    td=8;
                    u = trnd(td,T,N);
                    y = zeros(T, N );
                    X = zeros(T, N, p);

                    for t=1:T
                        if t <= T1
                            a0 = a0A;              % before break, a0A
                            beta0 = beta0A;
                        else
                            a0 = a0B;              % after break, a0B
                            beta0 = beta0B;
                        end

                        for i=1:N
                             % generate xit, uit, yit
                            X(t,i,1:p) = mu0(i,1)+ random('Normal', 0, sigma_x, 1, p);
                            y(t,i) =  mu0(i,1)+reshape( X(t,i,1:p), [1 p]) * beta0(i, 1:p)' + u(t,i);
 
                        end
                    end
                    %% reshape
                    y = reshape(y, N*T, 1);
                    X = reshape(X, N*T, p);
                    u = reshape(u, N*T, 1);

                    clear beta0 beta0A beta0B u mu0
                    %% Save result
                    DataFolder = fullfile( sprintf(append('DATA/DATA_FE_Tdist_Type', num2str( Type ),'_BF',num2str(Bfrac),'_N',num2str( N ),'T', num2str( T )  )) );
                    if ~exist(DataFolder, 'dir')
                        mkdir(DataFolder)
                    end

                    name = [DataFolder,'/DATA_',num2str(ii)];

                    save(name)


                    save ii ii N T Type Bfrac
                    clear 
                    load ii ii N T Type Bfrac


                    toc;
                end

            end
        end
    end
end




