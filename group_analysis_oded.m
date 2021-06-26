
%% load into the variable name "data"
DATA_DIMS = 3;
data = load(sprintf('grps_data_%dD.mat',DATA_DIMS)).('grps_data');
%% 
CONDITION_NAMES = {'HC','FC','acq'};
% group numbers across conditions:

group_nums_hc = [];
n_neurons_in_grp_hc = [];
neurons_locations_hc = {};
group_nums_fc = [];
n_neurons_in_grp_fc = [];
neurons_locations_fc = {};
group_nums_aq = [];
n_neurons_in_grp_aq = [];
neurons_locations_aq = {};
for i =1:length(data)
    if contains(data(i).cond_name,'HC')
        group_nums_hc = [group_nums_hc data(i).num_of_groups];
        n_neurons_in_grp_hc = [n_neurons_in_grp_hc data(i).grp_sizes];
        neurons_locations_hc{size(neurons_locations_hc,2)+1} = data(i).cell_coord;
    elseif contains(data(i).cond_name,'FC')
        group_nums_fc = [group_nums_fc data(i).num_of_groups];
        n_neurons_in_grp_fc = [n_neurons_in_grp_fc data(i).grp_sizes];
        neurons_locations_fc{size(neurons_locations_fc,2)+1} = data(i).cell_coord;
    else
        group_nums_aq = [group_nums_aq data(i).num_of_groups];
        n_neurons_in_grp_aq = [n_neurons_in_grp_aq data(i).grp_sizes];
        neurons_locations_aq{size(neurons_locations_aq,2)+1} = data(i).cell_coord;
    end
end
%%
% test for differences in group numbers
ttest2(group_nums_hc,group_nums_fc)
ttest2(group_nums_hc,group_nums_aq)
ttest2(group_nums_aq,group_nums_fc)


% test for differences in group sizes
ttest2(n_neurons_in_grp_hc,n_neurons_in_grp_fc)
ttest2(n_neurons_in_grp_hc,n_neurons_in_grp_aq)
ttest2(n_neurons_in_grp_aq,n_neurons_in_grp_fc)


% test for mutual connectivity
% The measure is for every neuron within a group,
% and is the number of neurons in the group that are closer than the
% threshold to that neuron, divided by the overall number of neurons in the
% group minus 1 (to exclude the neuron itself).
% This measure can than be used to give the average "clusteriness" of a
% group, as well as the standard deviation of that group in this measure
threshold=30*1.5;
[hc_conn_score,hc_conn_std]=get_condition_group_connectivity_metrics(threshold,neurons_locations_hc,DATA_DIMS);
[fc_conn_score,fc_conn_std]=get_condition_group_connectivity_metrics(threshold,neurons_locations_fc,DATA_DIMS);
[aq_conn_score,aq_conn_std]=get_condition_group_connectivity_metrics(threshold,neurons_locations_aq,DATA_DIMS);


ttest2(hc_conn_score,fc_conn_score)
ttest2(hc_conn_score,aq_conn_score)
ttest2(fc_conn_score,aq_conn_score)

ttest2(hc_conn_std,fc_conn_std)
ttest2(hc_conn_std,aq_conn_std)
ttest2(fc_conn_std,aq_conn_std)
