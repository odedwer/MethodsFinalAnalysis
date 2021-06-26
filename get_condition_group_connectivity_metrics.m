function [scores,score_stds] = get_condition_group_connectivity_metrics(threshold, cond_cell_arr,dims)
%GET_CONDITION_GROUP_CONNECTIVITY_METRICS
% This function calculates the mean and STD of connectivity scores of each group in the
% given condition. Read about the measure in get_connectivity_score
%
% threshold - the threshold to calculate the score by
% cond_cell_arr - cell array of groups. Each cell contains a cell array of
%       arrays with grouped neurons locations. This means cond_cell_arr{2}{3} is
%       a (N x dims) matrix  containing the coordinates in dims dimensions  
%       of the N neurons group 3 in  subject 2
% dims - number of dimensions to use, 1-3 inclusive
num_elements = 0;
for subj=1:size(cond_cell_arr,2)
    num_elements=num_elements+size(cond_cell_arr{subj},1);
end
curr_idx = 0;
scores = zeros(1,num_elements);
score_stds = zeros(1,num_elements);
for subj=1:size(cond_cell_arr,2)
    for grp=1:size(cond_cell_arr{subj},1)
        curr_idx=curr_idx+1;
        switch dims
            case 1
                cur_scores = get_connectivity_scores(threshold,cond_cell_arr{subj}{grp}(:,1));
            case 2
                cur_scores = get_connectivity_scores(threshold,cond_cell_arr{subj}{grp}(:,1),...
                    cond_cell_arr{subj}{grp}(:,2));
            case 3
                cur_scores = get_connectivity_scores(threshold,cond_cell_arr{subj}{grp}(:,1),...
                    cond_cell_arr{subj}{grp}(:,2),cond_cell_arr{subj}{grp}(:,3));
            otherwise
                error(sprintf('dims must be 1,2, or 3, but is %d',dims))
        end
        scores(1,curr_idx) =mean(cur_scores);
        score_stds(1,curr_idx) = std(cur_scores);
    end
end

end

