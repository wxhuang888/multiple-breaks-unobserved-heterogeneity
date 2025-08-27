% This code implements the LumsdaineOkuiWang2023 LSGB estimation procedure 
% for DGPG 5

warning off
clc
clear all

% addpath(genpath('reference'));
addpath('reference/LumsdaineOkuiWang2023replicationcode/replication')

%% setting: specification
for M    = [2]                   % # breaks: 1,2
for N    = [100]           % 200, 400
for T    = [100]           % 200, 400
sigma_e = 1;

%% folders(REVISE YOUR FOLDERS HERE)
% folder to store data
DataFolder = '';
% DataFolder = fullfile( sprintf(append('../2_GroupPanel/DATA/DATA_FE_Type',num2str( Type ), '_N',num2str( N ),'T', num2str( T )  )) );
% DataFolder = fullfile( sprintf(append('DATA/DATA_FE_Tdist_Multi_N',num2str( N ),'T', num2str( T )  )) );

% folder to store results
ResultFolder = '';
% ResultFolder = fullfile( sprintf(append('RESULTCOMP_LSGB/RESULT_G_FE_Type',num2str( Type ), ...
%                                            '/RESULT_M',num2str( M ),'_N',num2str( N ),'T', num2str( T ), '_s2n', num2str(sigma_e)  )) );
% ResultFolder = fullfile( sprintf(append('RESULTCOMP_LSGB/RESULT_G_FE_Multi', ...
%                                            '/RESULT_M',num2str( M ),'_N',num2str( N ),'T', num2str( T ), '_s2n', num2str(sigma_e)  )) );
if ~exist(ResultFolder, 'dir')
    mkdir(ResultFolder)
end

%% run 500 replications
for i=1
    i
    FILENAME    = [DataFolder,'/DATA_', num2str(i),'.mat'];
    NEWFILENAME = [ResultFolder,'/RESULT_', num2str(i),'.mat' ];

    %% estimation of T0
    if ~isfile(NEWFILENAME)
    tic;  % 20s (N50T50)
    try
        EST_comp_lsgb_seqM2(FILENAME, NEWFILENAME); %Lumsdaine, Okui, and Wang (2023)
    catch
        sprintf('%d has THE error',i);
    end
    toc;
    end
    
end


end
end
end



% end


