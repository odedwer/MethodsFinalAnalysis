%% Create structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts takes the 2D / 3D LabeledPoints mat files and generates a
% structure with the following components for each condition (=separate
% original mat files: 12 for 3D, 30 for 2D):
% cond_name: The name of the original mat file (e.g., "FC1_000_100_R_hemi_neuro_spots.mat").
% num_of_groups: The number of groups yielded by
% grp_shapes: the dbscan algorithm.
% grp_idx: a vector holding all of the group numers assigned by dbscan.
% cell_idx: the indices (taken from the LabeledPoints file) of cells
% clustered to each group (cell array).
% cell_coord: The [x,y,z] coordinates of the corresponding clustered cells (cell array).
% grp_shapes: The alphashapes (see https://www.mathworks.com/help/matlab/ref/alphashape.html) 
% of each group (cell array).
% grp_total_lngths: The area / volume of each group's alphashape (2D or 3D
% data). Note that for 3D, at least 4 cells are necessary to calculate the volume, 
% whereas for 2D at least 3 cels are necessary. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
choose_data='3D';
LabeledPoints = load(sprintf('LabeledPoints_%s.mat',choose_data)).(sprintf('LabeledPoints%s',choose_data));
grps_data=[];
cond_names=unique(LabeledPoints(:,5));
num_of_cond=length(cond_names);
for i=1:num_of_cond
    neuron_coord=[];
    shp=[];
    total_lngths=[];
    cond_handle=cond_names(i);
    cond_ind=find(ismember(LabeledPoints(:,5),cond_handle)==1);
    grp_idx=unique(double(LabeledPoints(cond_ind,4)));
    if isempty(find(ismember(grp_idx,-1)==1))==0
        grp_idx(find(ismember(grp_idx,-1)==1))=[];
    else
    end
    num_of_grps=length(grp_idx);
    neuron_ind=cell(num_of_grps,1);
    for j=1:num_of_grps
        neuron_ind{j}=cond_ind(find(double(LabeledPoints(cond_ind,4))==grp_idx(j)));
        if strcmp(choose_data,'3D')
            neuron_coord{j,1}=str2double(LabeledPoints(neuron_ind{j},1:3));
            for k=1:length(neuron_coord)
                shp{j,1}=alphaShape(neuron_coord{j,1},Inf);
                total_lngths{j,1}=volume(shp{j,1});
            end
        elseif strcmp(choose_data,'2D')
            neuron_coord{j,1}=str2double(LabeledPoints(neuron_ind{j},1:2));
            for k=1:length(neuron_coord)
                shp{j,1}=alphaShape(neuron_coord{j,1},Inf);
                total_lngths{j,1}=area(shp{j,1});
            end
        end
    end    
    grps_data(i).cond_name=cond_handle;
    grps_data(i).num_of_groups=num_of_grps;
    grps_data(i).grp_idx=grp_idx;
    grps_data(i).cell_idx=neuron_ind;
    grps_data(i).grp_sizes = zeros(1,num_of_grps);
    for j=1:num_of_grps
        grps_data(i).grp_sizes(j)=length(neuron_ind{j});
    end
    grps_data(i).cell_coord=neuron_coord;
    grps_data(i).grp_shapes=shp;
    grps_data(i).grp_total_lngths=total_lngths;   
end
  
 save(sprintf('grps_data_%s.mat',choose_data),'grps_data')

%% Visualize clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script plots the cells (as spheres using scatter or scatter3)
% and alphashapes of each group for a given condition (number within the
% grps_data.mat specified by the variable cond_num).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
choose_data='3D';
load(sprintf('LabeledPoints_%s.mat',choose_data))
load(sprintf('grps_data_%s.mat',choose_data))
cond_num=11;
num_gr=grps_data(cond_num).num_of_groups;
for i=1:num_gr
    coordinates=cell2mat(grps_data(cond_num).cell_coord(i));
    shp=grps_data(cond_num).grp_shapes(i);
    if strcmp(choose_data,'3D')
        scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3))
    elseif strcmp(choose_data,'2D')
        scatter(coordinates(:,1),coordinates(:,2))
    end
    hold on
    plot(shp{1},'FaceColor',i*[1/num_gr,1/num_gr,1/num_gr],'FaceAlpha',0.1)
    hold on
end
axis equal

%% Ellipsoid fitting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script fits groups with at least 9 cells to ellipsoids. See script of ellipsoid_fit_new.m 
% for more explanations.
% con_num: Specify condition number
% gr_num: Specify group number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% con_num=;
% gr_num=;
% [ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new(cell2mat(grps_data(con_num).cell_coord(gr_num)));