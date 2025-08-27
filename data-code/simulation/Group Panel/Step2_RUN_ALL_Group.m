
% This code implements the estimation procedure 
% for DGPG 1--4


warning off
clc
clear all

pathhtc = '';
%% setting: specification
for Type    = [1]             % DGPG: 1,2,3,4
    for N    = [50 ]          % 50, 100
        for T    = [100]          % 50, 100
            %% setting for estimation: tunings
            M=1;
            
            % r1      = ⌊cT^{c1}⌋
            est.c  = 1;
            est.c1 = 0.6;   
            
            % lambda0 = const * N^(est.cN) * T^(est.cT)
            est.clambda0 = 1;
            est.cN = 1;
            est.cT = 0.15;

            %% folders(REVISE YOUR FOLDERS HERE)
            % folder to store data
             DataFolder = fullfile( sprintf(append(pathhtc, 'DATA/DATA_FE_Tdist_Type', num2str( Type ),'_N',num2str( N ),'T', num2str( T )  )) );

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

                %% estimation of T0
                tic;
               
                try
                    EST_Hetero_Inst_rolling(FILENAME, NEWFILENAME, est);
                    EST_Group_Refine_FERobust_NEW(FILENAME, NEWFILENAME)
                catch
                    sprintf('%d has THE error',i);
                end
                
                toc;

            end


        end
    end
end









