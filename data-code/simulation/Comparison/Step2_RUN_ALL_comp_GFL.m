% This code implements the Qian and Su's (2016) GFL estimation procedure 
% this needs to be run on Windows

warning off
clc
clear all

addpath(genpath('reference/GFL1.0'));
%% setting: specification
NF=0.5;
for DGPnum=[1 ]                % 
for M    = [1 ]                % # breaks: 0,1,2
for N    = [50 100 ]           % 50, 100
for T    = [50 100 ]           % 50, 100
sigma_e = 1;

%% folders(REVISE YOUR FOLDERS HERE)
% folder to store data
DataFolder = '';   
% DataFolder = fullfile( sprintf(append('DGPs_Hetero_DATA_FE_U/DATA/DATA', num2str( DGPnum ), '_U_M',num2str( M ),'_NF',num2str( NF ),'_N',num2str( N ),'_T', num2str( T )  )) );

% folder to store results
ResultFolder = '';  
% ResultFolder = fullfile( sprintf(append('RESULTCOMP_GFL/RESULT_D', num2str( DGPnum ), ...
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
        EST_comp_gfl(FILENAME, NEWFILENAME);    
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



% end


