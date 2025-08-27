
% This code implements the estimation procedure generate 


warning off
clc
clear all

pathhtc = '';
%% setting: specification
for const =[1]  % constant for tuning parameter lambda0 = const * N^(est.cN) * T^(est.cT)
    for c = [1] % constant for tuning parameter r1      = ⌊cT^{c1}⌋
        for DGPnum=[1 ]       % DGPI: 1,2,3,4
            for M    = [2 ]   % # breaks
                for N    = [50 ]          
                    for T    = [100 ]         
                        %% setting: for DGPs
                        NF=0.5;
                        jmin=1;
                        jmax=2;
                        
                        %% setting: tunings for estimation
                        % r1      = ⌊cT^{c1}⌋
                        est.c  = c;
                        est.c1 = 0.6;   

                        % lambda0 = const * N^(est.cN) * T^(est.cT)
                        est.clambda0 = const;
                        est.cN = 1;
                        est.cT = 0.15;

                        %% folders(REVISE YOUR FOLDERS HERE)
                        % folder to store data
                        DataFolder = fullfile( sprintf(append(pathhtc, 'DATA/DATA',num2str( DGPnum ),'_U_M',num2str( M ),'_NF',num2str( NF ),'_N',num2str( N ),'_T', num2str( T ) )) );

                        % folder to store results
                        ResultFolder = fullfile( sprintf(append(pathhtc, 'RESULT/RESULT_U_DATA',num2str( DGPnum ),'_NF',num2str( NF ), '_c_', num2str( est.c ), '_c1_', num2str( est.c1 ), '_clambda_', num2str( est.clambda0 ), '_cN_', num2str( est.cN ), '_cT_', num2str( est.cT ), ...
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
                            % if ~isfile(NEWFILENAME)
                                try
                                    EST_Hetero_Inst_rolling(FILENAME, NEWFILENAME, est);
                                catch
                                    sprintf('%d has THE error',i);
                                end
                            % end
                            toc;

                        end


                    end
                end
            end
        end
    end
end









