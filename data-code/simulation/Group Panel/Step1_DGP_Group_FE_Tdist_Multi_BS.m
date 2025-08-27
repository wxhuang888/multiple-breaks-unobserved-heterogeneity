
% This code generates Group Panel with M breaks in slope coefficient \beta
% for DGPG 5 
% for different distances between two breaks (T_2^0 - T_1^0)/T


warning off
clc
clear all

%% Replicate 500
for BS=[0.1:0.05:0.5]      % (T_2^0 - T_1^0)/T
for N = [50 100]
    for T = [50 100]
        for ii=1:1
            disp(ii)
            
            tic;
            %% setting

            RandomB1=randn(1);
            RandomB2=randn(1);
            T1 = floor(1/3 * T+RandomB1);      % break point 1
            T2 = floor(T1+BS*T);      % break point 2
            regime0 = [1;T1+1;T2+1;T+1];
            K = 2;    % # of groups
            p = 1;    % # of regressors

            % two breaks
            N1A = 0.3*N;
            N1B = N1A;
            N1C = 0.7*N;

            a0A = [-0.4;
                1.6];
            a0B = [0.4;
                0.8];
            a0C = [-0.4;
                1.6];


            % same as a0A, a0B, but for each series
            beta0A = [ones(N1A,1) * a0A(1,:);
                ones((N-N1A),1) * a0A(2,:)];

            beta0B = [ones(N1B,1) * a0B(1,:);
                ones((N-N1B),1) * a0B(2,:)];

            beta0C = [ones(N1C,1) * a0C(1,:);
                ones((N-N1C),1) * a0C(2,:)];

            % signal - noise
            sigma_x = 1;
            sigma_e = 1;

            %% DGP

            mu0=randn(N,1);
            td=8;
            u = sigma_e * trnd(td,T,N);
            y = zeros(T, N );
            X = zeros(T, N, p);

            for t=1:T
                if t <= T1
                    a0 = a0A;              % before break 1, a0A
                    beta0 = beta0A;
                elseif t <= T2
                    a0 = a0B;              % after break 1, a0B
                    beta0 = beta0B;
                else
                    a0 = a0C;              % after break 2, a0C
                    beta0 = beta0C;
                end

                for i=1:N

                    % generate xit, uit, yit
                    X(t,i,1:p) =  mu0(i,1)+random('Normal', 0, sigma_x, 1, p);
                    y(t,i) = mu0(i,1) + reshape( X(t,i,1:p), [1 p]) * beta0(i, 1:p)' + u(t,i);

                end
            end

            % reshape
            y = reshape(y, N*T, 1);
            X = reshape(X, N*T, p);
            u = reshape(u, N*T, 1);

            clear beta0 beta0A beta0B beta0C u mu0
            %% Save result
            DataFolder = fullfile( sprintf(append('DATA/DATA_FE_Tdist_Multi', '_BF', num2str(BS), '_N',num2str( N ),'T', num2str( T )  )) );
            if ~exist(DataFolder, 'dir')
                mkdir(DataFolder)
            end

            name = [DataFolder,'/DATA_',num2str(ii)];

            save(name)


            save ii ii N T BS
            clear all
            load ii ii N T BS

            toc;
        end
    end
end
end




