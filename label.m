%% introduction
% run this script in a folder containing the data
% (it's okay if the *.mat files are in subfolders)
% it will create a LabeledPoints.mat file containing a matrix
% first 3 columnts are x,y,z coordinates 
% 4'th column is the label if the group
% 5'th column is the file name (to know the source of the date) 

%% global arguments
myFolder = '.';
inflation_threshold = 1.5;
epsilon_2D = 30;
epsilon_3D = epsilon_2D*inflation_threshold;
minpts = 2;
filePattern = fullfile(myFolder, '**/*_*.mat');
MatFiles = dir(filePattern);
LabeledPoints2D = [];
LabeledPoints3D = [];
%%
for mat = 1 : length(MatFiles)
    % read matrix
    M = load([MatFiles(mat).folder '/' MatFiles(mat).name],'*');
    M = cell2mat(struct2cell(M));
    NumOfNeurons = size(M,1);
    FileLabel = strings(NumOfNeurons,1);
    FileLabel(1:NumOfNeurons) = string(MatFiles(mat).name);
    
    % dbscan
    if (contains(MatFiles(mat).folder,'2D'))
        idx = dbscan(M,epsilon_2D, minpts);
        M = [M idx FileLabel];
        LabeledPoints2D = [LabeledPoints2D;M];
    elseif(contains(MatFiles(mat).folder,'3D'))
        idx = dbscan(M,epsilon_3D, minpts);
        M = [M idx FileLabel];
        LabeledPoints3D = [LabeledPoints3D;M];
    end
end
%% save
save('LabeledPoints_2D.mat','LabeledPoints2D');
save('LabeledPoints_3D.mat','LabeledPoints3D');