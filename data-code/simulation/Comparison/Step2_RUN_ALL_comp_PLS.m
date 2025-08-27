% This code implements the Qian and Su's (2016) PLS estimation procedure 

warning off
clc
clear all

addpath(genpath('reference'));
pathhtc='';
%% setting: specification
for DGPnum = [1];       
for M    = [1 2]              % # breaks: 1,2
for N    = [50 100]           % 50, 100
for T    = [50 100]           % 50, 100
sigma_e = 1;

%% folders(REVISE YOUR FOLDERS HERE)
% folder to store data
DataFolder = '';
% DataFolder = fullfile( sprintf(append(pathhtc, 'DATA/DATA',num2str( DGPnum ),'_U_RB_M',num2str( M ),'_N',num2str( N ),'T', num2str( T ))) );


% folder to store results
ResultFolder = '';
% ResultFolder = fullfile( sprintf(append('RESULTCOMP_PLS_U/RESULT_D', num2str( DGPnum ), ...
%                                            '/RESULT_M',num2str( M ),'_N',num2str( N ),'T', num2str( T )  )) );
if ~exist(ResultFolder, 'dir')
    mkdir(ResultFolder)
end

%% run 500 replications
for i=1:1
    i
    FILENAME    = [DataFolder,'DATA_', num2str(i),'.mat'];
    NEWFILENAME = [ResultFolder,'RESULT_', num2str(i),'.mat' ];

    %% estimation of T0
    % tic;
    % if ~isfile(NEWFILENAME)
    try
        EST_comp_pls(FILENAME, NEWFILENAME);
    catch
        sprintf('%d has THE error',i);
    end
    % end
    % toc;
    
end


end
end
end
end


