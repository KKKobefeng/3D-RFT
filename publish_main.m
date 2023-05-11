clear
close all

%% Publish RFT Functions
% List the files you want to include in the published output, except main.m
file_list = {
    'src/rotateTriangulationX.m',
    'src/moveTriangulationZ.m',
    'src/getStlData.m',
    'src/calcVelocity.m',
    'src/checkConditions.m',
    'src/findLocalFrame.m',
    'src/findAngles.m',
    'src/findFit.m',
    'src/findAlpha.m',
    'src/getForces.m',
    'src/getTorques.m'
};

% Create a temporary file to store the concatenated contents
temp_file = 'appendix_RFTfunc.m';
fid = fopen(temp_file, 'w');

% Concatenate the contents of all the files in file_list
for i = 1:length(file_list)
    file_contents = fileread(file_list{i});
    fprintf(fid, '%%%% Function: %s\n\n%s\n\n', file_list{i}, file_contents);
end

% Close the temporary file
fclose(fid);

options_others = struct('format', 'pdf', 'outputDir', 'output', 'evalCode', false);

% Publish the temporary file without code execution
publish(temp_file, options_others);

% Delete the temporary file
delete(temp_file);


%% Publish Main Function
% Set the options for the publish function
options_main = struct('format', 'pdf', 'outputDir', 'output', 'evalCode', true, ...
                      'figureSnapMethod', 'entireGUIWindow', 'imageFormat', 'bmp', 'maxHeight', 300);

% Publish the main.m file with code execution enabled
publish('main.m', options_main);

% figure_files = dir(fullfile('output', '*.png'));
% for i = 1:length(figure_files)
%     delete(fullfile('output', figure_files(i).name));
% end
