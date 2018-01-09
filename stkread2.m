function [ stack ] = stkread2( filename )
%STKREAD2 reads the entire tiff stack
%   [ stack ] = stkread2( filename )

if nargin < 1
    [filename, filepath] = uigetfile('*.tif');
    filename = fullfile(filepath, filename);
end


% Obtain the number of frames
n_frames = length(imfinfo(filename));

% Read first frame
stack = imread(filename, 1);

% Expand the stack
stack = repmat(stack, [ 1 1 n_frames]);

% Load the stack
for i = 2 : n_frames
    
    % Load the frames
    stack(:,:,i) = imread(filename, i);
    
end

end

