% move_ct_files.m
% This script takes all of the series in a large CT scan exam and sorts
% them by series (scan number).

clear all
close all

% Path to folder with CT Data
cd('E:\Experiment Data...')% Series number to put in folders
% Scan numbers
series = [1:3];

% Number of images per series
images_per_series = 88;
% Prefix for save folder names, it is easiest to put nothing
file_prefix = '';

% Look at image strings and determine position of space between series
% number and image number
series_number_space = 57;
% Loop through each series
for i = 1:length(series)
    % Display update on sorting
    disp(['Currently sorting files in series ', num2str(series(i))])
    % Make directory for series(i)
    mkdir([file_prefix, num2str(series(i))])
    dir_list = dir;
    % Account for fixed string length of series number in CT filename
    if series(i) < 100
        ct_series_number = ['00',num2str(series(i)),'_'];
    elseif series(i) > 99
        ct_series_number = ['0',num2str(series(i)),'_'];
    end
    
    % Search through files in directory for images with series number i
    for j =1:length(dir_list)
        str_match_test = strfind(dir_list(j).name, ct_series_number);
        % If there is a match
        if length(str_match_test) >= 1 
            % If that match is with the series number, NOT image number
            % then move the file to the folder created for that series
            if min(str_match_test)< series_number_space
                movefile([cd, '\',dir_list(j).name],[cd,'\', file_prefix, num2str(series(i))])
            end
        end
    end 
    
    % Verify that there are required number of files in folder
    dir_check = dir([[cd,'\', file_prefix, num2str(series(i))], '\*.v2']);
    if images_per_series ~= length(dir_check)
        disp(['WARNING: There are missing images in folder ',  ...
            file_prefix, num2str(series(i))])
    end
end
