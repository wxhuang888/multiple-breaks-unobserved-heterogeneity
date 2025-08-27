% This code implements the BFK 2016 estimation procedure 

warning off
clc
clear all

addpath(genpath('reference'));
%% setting: specification
for DGPnum = 1                 % DGP1~4
    for M    = [1 2]                  % # breaks: 1,2
        for N    = [50 100]           % 200, 400
            for T    = [50 100]           % 200, 400
                sigma_e = 1;

                %% folders (REVISE YOUR FOLDERS HERE)
                % folder to store data
                DataFolder = '';   
                
                % folder to store results
                ResultFolder = '';  
                % if ~exist(ResultFolder, 'dir')
                %     mkdir(ResultFolder)
                % end

                %% run 500 replications
                for i=1:1
                    i
                    FILENAME    = [DataFolder,'/DATA_', num2str(i),'.mat'];
                    NEWFILENAME = [ResultFolder,'/RESULT_', num2str(i),'.mat' ];

                    %% estimation of T0
                    tic;
                    try
                        EST_comp_bfk_Mknown(FILENAME, NEWFILENAME, M);
                    catch
                        sprintf('%d has THE error',i);
                    end
                    toc;

                end


            end
        end
    end
end





