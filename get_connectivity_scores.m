function [scores] = get_connectivity_scores(threshold,x,y,z)
%GET_CONNECTIVITY_SCORES
% this function calculates the connectivity score of the given points
% threshold - the threshold under which we consider points to be connected
% x - the x values of the points, a row vector
% y - Optional, the y values of the points, a row vector
% z - Optional, the z values of the points, a row vector
dist_mat = (x-x').^2;
if nargin==3
    dist_mat = dist_mat+(y-y').^2;
end
if nargin==4
    dist_mat = dist_mat+(z-z').^2;
end
dist_mat = sqrt(dist_mat);
dist_mat(logical(eye(size(dist_mat,1))))=nan;
scores =  (nansum(dist_mat<=threshold,1))/(numel(x)-1);
end

