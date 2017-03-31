function [ optimalthreshold, areamin ] = thresholdgradient( inputimg, threshold_range, increment )
%THRESHOLDGRADIENT tries a range of thresholds to determine which one gives
%the most consistent ROI size (std/mean). It also returns the minimal area
%size.
%   [ optimalthreshold, areamin ] = thresholdgradient( inputimg, threshold_range )


if nargin < 3
    % Default increment
    increment = 10;
    
    if nargin < 2
        % Default range
        threshold_range = [10, 100];
    end
end

% Initiate thresholds to try
n_tries = abs(diff(threshold_range)/increment + 1);

% Thresholds to try
thresholds = threshold_range(1) + (0 : (n_tries-1)) * increment;
N_rois = zeros(n_tries, 1);
ROIsizestd = zeros(n_tries, 1);

% Min areas
area_min_vec = zeros(n_tries, 1);

for i = 1 : n_tries
    
    % Generate the number of ROIs
    [labeledimg, N_rois(i)] = bwlabel(inputimg >= thresholds(i));
    
    % Get areas
    roiproperties = regionprops(labeledimg);
    Areas = cell2mat({roiproperties.Area});
    
    if ~isempty(Areas)
        % Calcualte area std/mean and minimal values
        ROIsizestd(i) = std(Areas)/mean(Areas);
        area_min_vec(i) = min(Areas);
    end
end

scatter(thresholds, ROIsizestd, 'o')
ylabel('Area Std/Mean')
xlabel('Threshold')

% Find first threshold value where std = mean
best_ind = find(ROIsizestd < 1, 1, 'first');

% Outputs
optimalthreshold = thresholds(best_ind);
areamin = area_min_vec(best_ind);

% figure
% imshow(inputimg >= optimalthreshold)

end

