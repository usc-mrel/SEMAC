function demo_convert_twix_to_ismrmrd(data_directory)

setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% demo_convert_twix_to_ismrmrd.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/23/2021, Last modified: 05/23/2021

%% Clean slate
% close all; clear; clc;

%% Define data directory
% data_directory ='F:\USC\MREL\LowField\LungImaging\T2measurement\Data\VOL139_BL\RL\RAW';
cd(data_directory)

%% Make folder
h5_folder = [data_directory,'/h5'];
noise_folder = [data_directory,'/noise'];
if ~exist(h5_folder)
    mkdir h5
end
if ~exist(noise_folder)
    mkdir noise
end

%% Get directory info
dir_info = dir(fullfile(data_directory, '*.dat'));

%% Convert TWIX to ISMRMRD format
nr_dat_files = length(dir_info);

start_time = tic;
for idx = 1:nr_dat_files
    dat_filename = dir_info(idx).name;
    dot_loc = strfind(dat_filename, '.');
    
    %----------------------------------------------------------------------
    % noise measurements
    % e.g.) Protocol name [1]: AdjQuietCoilSens
    %----------------------------------------------------------------------
    linux_command1 = sprintf('siemens_to_ismrmrd -f %s -z 1 -o noise_%s.h5', dat_filename, dat_filename(1:dot_loc-1));
    tic; fprintf('(%2d/%2d): Running %s... ', idx, nr_dat_files, linux_command1);
    [status1,result1] = system(linux_command1);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));    
    
    %----------------------------------------------------------------------
    % k-space measurements
    % e.g.) Protocol name [2]: tse_ES_RS10_n15000
    %----------------------------------------------------------------------
    linux_command2 = sprintf('siemens_to_ismrmrd -f %s -z 2 -o %s.h5', dat_filename, dat_filename(1:dot_loc-1));
    tic; fprintf('(%2d/%2d): Running %s... ', idx, nr_dat_files, linux_command2);
    [status2,result2] = system(linux_command2);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc, toc(start_time));
end

movefile noise_* noise
movefile *.h5 h5