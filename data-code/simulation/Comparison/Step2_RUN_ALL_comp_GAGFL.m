% This code implements the OkuiWang2021 GAGFL estimation procedure 

warning off
clc
clear all

% addpath(genpath('../reference'));
addpath('reference/OkuiWang2021replicationcode/replication files/simulation_replication')

%% setting: specification
Type = 1;                     % 1,2,3,4,multi
for M    = [1 ]               % # breaks: 1,2
for N    = [50 100]           % 200, 400
for T    = [50 100]           % 200, 400
sigma_e = 1;
G = 2;          % # groups
%% folders
% folder to store data
DataFolder = '';
% DataFolder = fullfile( sprintf(append('DATA/DATA_FE_Tdist_Type',num2str( Type ), '_N',num2str( N ),'T', num2str( T )  )) );

% folder to store results
ResultFolder = '';
% ResultFolder = fullfile( sprintf(append('RESULTCOMP_GAGFL/RESULT_G_FE_Type',num2str( Type ), ...
%                                            '/RESULT_M',num2str( M ),'_N',num2str( N ),'T', num2str( T ), '_s2n', num2str(sigma_e)  )) );
if ~exist(ResultFolder, 'dir')
    mkdir(ResultFolder)
end

%% run 500 replications
for i=1:1
    i
    FILENAME    = [DataFolder,'/DATA_', num2str(i),'.mat'];
    NEWFILENAME = [ResultFolder,'/RESULT_', num2str(i),'.mat' ];

    %% estimation of T0
    if ~isfile(NEWFILENAME)
    % tic;  % 2~3min (N50T50)
    try
        EST_comp_gagfl(FILENAME, NEWFILENAME, G);   % Okui and Wang (2021)
    catch
        sprintf('%d has THE error',i);
    end
    % toc;
    end
    
end


end
end
end



% end


