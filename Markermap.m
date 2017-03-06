function [  ] = Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh, M2thresh, NROIs_master)
%MARKERMAP color-codes the dapi segmentation based on marker labelings
%   Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh, M2thresh, NROIs_master)

nslices = size(Dapireg,3);

hwait = waitbar(0, 'Generating map');

% Generate RGB
% RGB2show = repmat(Dapireg(:,:,ii), [1 1 3]);
% RGB2show(:,:,1) = 0;
% RGB2show(:,:,2) = 0;

M1temp = zeros(size(Dapireg(:,:,1)));
M2temp = zeros(size(Dapireg(:,:,1)));

% Make a map across slices
M1map = repmat(M1temp, [1 1 nslices]);
M2map = repmat(M2temp, [1 1 nslices]);

for ii = 1 : nslices
    waitbar(ii/nslices)


    % Determine positivity

    M1positive =...
        M1pix_cell{ii} >= M1thresh(ii);
    M2positive =...
        M2pix_cell{ii} >= M2thresh(ii);
 

    % Color-code the positivity
    for i = 1 : NROIs_master(ii)

        % Marker 1
        M1temp = M1temp + ...
            (Dapireg(:,:,ii) == i) * double(M1positive(i));

        % Marker 2
        M2temp = M2temp + ...
            (Dapireg(:,:,ii) == i) * double(M2positive(i));

    end
    
    % Store maps
    M1map(:,:,ii) = M1temp;
    M2map(:,:,ii) = M2temp;
end

close(hwait)

Dapireg_max = max(Dapireg,[],3);
RGB2show = repmat(Dapireg_max , [1 1 3]);
RGB2show(:,:,1) = max(M1map,[],3);
RGB2show(:,:,2) = max(M2map,[],3);

imshow(RGB2show)

end

