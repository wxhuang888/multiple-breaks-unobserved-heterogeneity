
% This code obtains the oracle group-specific coefficient
% for DGPG 1--4


warning off
clc
clear all

pathhtc = '';
%% setting: specification
for Type =[1 ]              % 1,2,3,4
    for M    = [1]           % # breaks: 0,1,2
        for N    = [50 ]          % 200, 400
            for T    = [100]          % 200, 400
                sigma_e = 1;
                %% setting: estimation
                est.c  = 1;
                est.c1 = 0.6;   

                % lambda0 = const * N^(est.cN) * T^(est.cT)
                est.clambda0 = 1;
                est.cN = 1;
                est.cT = 0.15;

                %% folders(REVISE YOUR FOLDERS HERE)
                % folder to store data
                DataFolder = fullfile( sprintf(append('DATA/DATA_FE_Tdist_Type', num2str( Type ),'_N',num2str( N ),'T', num2str( T )  )) );

                % folder to store results
                ResultFolder = fullfile( sprintf(append(pathhtc, 'RESULT/RESULT_Type', num2str( Type ), '_c_', num2str( est.c ), '_c1_', num2str( est.c1 ), '_clambda_', num2str( est.clambda0 ), '_cN_', num2str( est.cN ), '_cT_', num2str( est.cT ), ...
                    '/RESULT_M',num2str( M ),'_N',num2str( N ),'T', num2str( T )  )) );
                if ~exist(ResultFolder, 'dir')
                    mkdir(ResultFolder)
                end

                
                %% run 500 replications
                for i=1:1
                    i
                    FILENAME    = [DataFolder,'/DATA_', num2str(i),'.mat'];
                    NEWFILENAME = [ResultFolder,'/RESULT_', num2str(i),'.mat' ];

                    %% estimation
                    tic;
                    try
                        EST_Oracle_Coeff(FILENAME, NEWFILENAME, Type);
                    catch
                        sprintf('%d has THE error',i);
                    end
                    toc;

                end


            end
        end
    end
end
