% The intensities in the first and the last slices of the stack would not be 
% quantified (they are zeros). You have to mannually delete them.

%% Parameters

Smooththresh = 0; % Determine if smooth threshold or not
Automarkerthresh = 0; % If automatedly determine marker threshold

firstslice2read = 1; % First slice to read
maxnslices = 22; % Only consider first X slides.
openclose = 5; % Size for the dirts
backgroundsize = 12; % Size to use for background subtraction
overlaptreshold = 0.5; % Fraction overlap to be considered the same cell
overlapretro = 3; % Use the past N slices for overlap considerations
Marker1 = 'LMO';
Marker2 = 'LHX1';

Commonfilapath = '/Users/Zhihua/Documents/Dropbox (HMS)/Confocal/20170210-SCN brain3-DAPI_LMO3_LMO3_LHX1/Brain 3 SCN quantify/Left';

%% Read files

% Grab file names
[Dapifile, filepath] = uigetfile(fullfile(Commonfilapath,'*.tif'), 'Dapi file:');
Marker1file = uigetfile(fullfile(filepath,'*.tif'), Marker1);
Marker2file = uigetfile(fullfile(filepath,'*.tif'), Marker2);

% Determine slice numbers
fileinfo = imfinfo(fullfile(filepath,Dapifile));
lastslice2read = min(maxnslices, size(fileinfo,1));
imgsize = [fileinfo(1).Width, fileinfo(1).Height];
nslices = lastslice2read - firstslice2read + 1;

% Prime the stacks
Dapistack = zeros(imgsize(2), imgsize(1), nslices);
Marker1stack = Dapistack;
Marker2stack = Dapistack;

% Prime the segmentation matrix
Dapireg = Dapistack;

% Read the files
for i = 1 : nslices
    Dapistack(:,:,i) = imread(fullfile(filepath, Dapifile), ...
        firstslice2read + i - 1);
    Marker1stack(:,:,i) = imread(fullfile(filepath, Marker1file), ...
        firstslice2read + i - 1);
    Marker2stack(:,:,i) = imread(fullfile(filepath, Marker2file), ...
        firstslice2read + i - 1);
end

% Replacing saturated dirts in LMO channel with background intensity 500

indices_LMO = Marker1stack >= 3000;
Marker1stack(indices_LMO) = 500;

% Code for testing
% figure; imshow(mat2gray(Marker1stack(:,:,4),[0,4095]))

%% Histogram scaling

%ZL Reduce the Heigth x Width x nSlice matrix to slice x nSlice.
Dapistack_2d = reshape(Dapistack,[],nslices);
Marker1stack_2d = reshape(Marker1stack,[],nslices);
Marker2stack_2d = reshape(Marker2stack,[],nslices);

%ZL remove the intensities of empty places (<400) and real signal (>4090)
indices_dapi = Dapistack_2d<480;
indices1 = Marker1stack_2d<300 | Marker1stack_2d>800;
indices2 = Marker2stack_2d<300 | Marker2stack_2d>600;
Dapistack_2d(indices_dapi) = NaN;
Marker1stack_2d(indices1) = NaN;
Marker2stack_2d(indices2) = NaN;

%ZL Get the most frequent intensity in each slices for background
%ZL Get the percentiles of the data set
%Mdapiscale = nanmedian(Dapistack_2d,1);
Mdapiscale = nanmedian(Dapistack_2d,1);
M1scale = nanmedian(Marker1stack_2d,1);
M2scale = nanmedian(Marker2stack_2d,1);

%ZL Rescale to the level of sclice with the max backgroud
[~,Mdapiscale_max] = max(Mdapiscale(:));
[~,M1scale_max] = max(M1scale(:));
[~,M2scale_max] = max(M2scale(:));
Mdapiscale = Mdapiscale/Mdapiscale(Mdapiscale_max);
M1scale = M1scale/M1scale(M1scale_max);
M2scale = M2scale/M2scale(M2scale_max);

% Plot scale for each channel
figure
plot(squeeze(Mdapiscale))
figure
plot(squeeze(M1scale))
figure
plot(squeeze(M2scale))

%Expand the dimensions
Mdapiscale = permute(Mdapiscale,[3,1,2]);
M1scale = permute(M1scale,[3,1,2]);
M2scale = permute(M2scale,[3,1,2]);

Mdapiscale_mat = repmat(Mdapiscale,[imgsize(2),imgsize(1),1]);
M1scale_mat = repmat(M1scale,[imgsize(2),imgsize(1),1]);
M2scale_mat = repmat(M2scale,[imgsize(2),imgsize(1),1]);

Dapistack_scaled = Dapistack ./ Mdapiscale_mat;
Marker1stack_scaled = Marker1stack ./ M1scale_mat;
Marker2stack_scaled = Marker2stack ./ M2scale_mat;

% figure;imshow(mat2gray(Dapistack_scaled(:,:,20),[0,4095]))
% figure;imshow(mat2gray(Marker1stack_scaled(:,:,20),[0,4095]))
% figure;imshow(mat2gray(Marker2stack_scaled(:,:,20),[0,4095]))
% figure; imshowpair(mat2gray(Dapistack_scaled(:,:,20)),mat2gray(Dapistack(:,:,20)),'montage')
% figure;imshow(mat2gray(Dapistack_scaled(:,:,20),[0,4095]))


%% Segmentate DAPI

Dapistack_scaled_nobg = zeros(size(Dapistack,1), size(Dapistack,2), nslices);
% Remove background in DAPI channel
for ii = 1 : nslices
    Dapibg = imopen(Dapistack_scaled(:,:,ii), strel('disk', backgroundsize));
    Dapistack_scaled_nobg(:,:,ii) = Dapistack_scaled(:,:,ii) - Dapibg;
end

% hwait = waitbar(0, 'Segmentating');
% dapithreshold = 250;
% areamin = 50;

for ii = 1 : nslices
    % waitbar(ii/nslices)
    
    % Open connections
    dapiopen = imopen(Dapistack_scaled_nobg(:,:,ii), strel('disk', openclose));
    
        % Calculate a threshol
        [thresholdvalue, areamin] = thresholdgradient(dapiopen, [200, 1000]);
    
    % Apply thresholds
    dapiroi = dapiopen >= thresholdvalue;
    
    % Watershed %ZL, change the Connectivity of Watershed (default 8)
    watershedded = watershed(-bwdist(~dapiroi));
    dapiroi_ws = dapiroi;
    dapiroi_ws(watershedded == 0) = 0;
    
    % Apply area threshold
    dapiroi_ws_ath = areathresh( dapiroi_ws, 0.5, areamin, 3, 8);
    
    % Store the registered image
    Dapireg(:,:,ii) = bwlabel(dapiroi_ws_ath);
    %     figure;imshowpair(Dapistack_scaled(:,:,ii),bwperim(Dapireg(:,:,ii)));
end

% close(hwait)

%% Remove overlaped ROIs across slices

% hwait = waitbar(0, 'Cleaning up overlap');

for   ii = 2 : nslices-1
    % Image to cleanup
    img2clean = Dapireg(:,:,ii);
    % Calculate how many cells are there in current slice
    n_areas = max(max(Dapireg(:,:,ii)));
    Dapipix = zeros(n_areas , 1);
    Dapi_current = Dapistack_scaled(:,:,ii);
    Dapi_previous=Dapistack_scaled(:,:,max(1,ii-1));
    Dapi_after=Dapistack_scaled(:,:,min(nslices,ii+1));
    
    % Examine all the cells one by one
    for i = 1 : n_areas
        roiarea = sum(sum(Dapireg(:,:,ii) == i));
        if mean(Dapi_current(Dapireg(:,:,ii) == i)) < mean(Dapi_previous(Dapireg(:,:,ii) == i)) ...
                || mean(Dapi_current(Dapireg(:,:,ii) == i)) < mean(Dapi_after(Dapireg(:,:,ii) == i)) ...
                || roiarea < 80
            img2clean(img2clean==i) = 0;
        end
    end
    % figure; imshow(img2clean);
    
    % Relabel and store
    Dapireg(:,:,ii) = bwlabel(img2clean>0);
   % figure;imshowpair(Dapistack_scaled(:,:,ii),bwperim(Dapireg(:,:,ii)));
end

% For testing, to examine overlaps between slices
for i = 2 : nslices-2
    figure;imshowpair(imadjust(Dapireg(:,:,i)),imadjust(Dapireg(:,:,i+1)));
end
%  figure;imshowpair(imadjust(Dapireg(:,:,1)),imadjust(Dapireg(:,:,2))); 
%  figure;imshowpair(Dapistack_scaled(:,:,12),bwperim(Dapireg(:,:,12)));
%  figure;imshowpair(Marker2stack_scaled(:,:,4),bwperim(Dapireg(:,:,4)));



%% Calculate pixel values for markers
% Prime pixel values for markers
% Prime a threshold vector
M1thresh = zeros(nslices,1);
M2thresh = zeros(nslices,1);
% Prime cells to contain pixel values
M1pix_cell = cell(nslices, 1);
M2pix_cell = cell(nslices, 1);

for ii = 2 : nslices-1
    % flatten backgrounds
    M1bg = imopen(Marker1stack_scaled(:,:,ii), strel('disk', backgroundsize));
    M2bg = imopen(Marker2stack_scaled(:,:,ii), strel('disk', backgroundsize));
    M1_nobg = Marker1stack_scaled(:,:,ii) - M1bg;
    M2_nobg = Marker2stack_scaled(:,:,ii) - M2bg;
    
    % Determine the number of ROIs
    n_areas = max(max(Dapireg(:,:,ii)));
    
    % Prime Pixel-value vector
    M1pix = zeros(n_areas , 1);
    M2pix = zeros(n_areas , 1);
    
    for i = 1 : n_areas
        % Calculate pixel values for both
        M1pix(i) = mean(M1_nobg(Dapireg(:,:,ii) == i));
        M2pix(i) = mean(M2_nobg(Dapireg(:,:,ii) == i));
    end
    
    % if ii == 1
    % Establish pixel intensity order
    [M1pix_sorted , M1order] = sort(M1pix, 1, 'descend');
    [M2pix_sorted , M2order] = sort(M2pix, 1, 'descend');
    
    %{  % Plot pixel intensities
    % Initialize parent
    %             fM1 = figure('name','Marker 1','Position',[650 100 600, 400]);
    %             plot(1:n_areas, M1pix_sorted,'-o');
    %             xlabel('Slices')
    %             grid on
    %
    %             fM2 = figure('name','Marker 2','Position',[1300 100 600, 400]);
    %             plot(1:n_areas, M2pix_sorted,'-o');
    %             xlabel('Slices')
    %             grid on
    
    % Use antomated or manual ways to determine threshold
    %     if Automarkerthresh == 1
    %         % mean + 2 std
    %         % M1thresh(ii) = mean(M1_nobg(:)) + 1.5 * std(M1_nobg(:));
    %         % M2thresh(ii) = mean(M1_nobg(:)) + 1.5 * std(M2_nobg(:));;
    %     else
    %}
    
    % Manual
    if M1thresh(ii) <= 0
        figurename = ['Slice ', num2str(ii), '/', num2str(nslices), ' - Marker 1'];
        % Manually find threshold if needed
        if max(max(Dapireg(:,:,ii))) > 0
            markerthreshold(Marker1stack_scaled(:,:,ii), Dapireg(:,:,ii),...
                M1order, figurename);
            uiwait()
            M1thresh(ii) = M1pix_sorted(min(lastpositive + 1, n_areas));
        end
    end
    
    if M2thresh(ii) <= 0
        figurename = ['Slice ', num2str(ii), ' - Marker 2'];
        
        % Manually find threshold if needed
        if max(max(Dapireg(:,:,ii))) > 0
            markerthreshold(Marker2stack_scaled(:,:,ii), Dapireg(:,:,ii),...
                M2order, figurename);
            uiwait()
            M2thresh(ii) = M2pix_sorted(min(lastpositive + 1, n_areas));
        end
    end
    
    %         close(fM1, fM2)
    
    % Load pixel values into cells
    M1pix_cell{ii} = M1pix;
    M2pix_cell{ii} = M2pix;
    
end

%% Smooth threshold
% if Smooththresh == 1 M1thresh_s = smooth(M1thresh); M2thresh_s = smooth(M2thresh); end

%% Consolidate data (Save data to master file)

% Total number of ROIs
NROIs_master = squeeze(max(max(Dapireg,[],1),[],2));

% Initiate matrix for total pixel values
M1pix_master = zeros(sum(NROIs_master,1),1);
M2pix_master = zeros(sum(NROIs_master,1),1);
M1positive_master = zeros(sum(NROIs_master,1),1);
M2positive_master = zeros(sum(NROIs_master,1),1);
Centroids_master = zeros(sum(NROIs_master,1),2);

% Book-keeping matrices
M1_cellind_master = zeros(sum(NROIs_master,1),1);
M1_slice_master = zeros(sum(NROIs_master,1),1);

% Loop through the slices
for ii = 1 : nslices
    % Determine start and end indices
    if ii == 1
        istart = 1;
    else
        istart = 1 + sum(NROIs_master(1 : (ii-1)));
    end
    
    iend = sum(NROIs_master(1:ii));
    
    if ii == 1 || ii == nslices
        M1pix_master(istart:iend) = 0;
        M2pix_master(istart:iend) = 0;
    else
        M1pix_master(istart:iend) = M1pix_cell{ii};
        M2pix_master(istart:iend) = M2pix_cell{ii};
    end
    
    % Determine positivity
    if ii == 1 || ii == nslices
        M1positive_master(istart:iend) = 0;
        M2positive_master(istart:iend) = 0;
    else
        M1positive_master(istart:iend) =...
            M1pix_cell{ii} > M1thresh(ii);
        M2positive_master(istart:iend) =...
            M2pix_cell{ii} > M2thresh(ii);
    end
    
%{  
       if Automarkerthresh < 1
%             if Smooththresh < 1
%                 M1positive_master(istart:iend) =...
%                     M1pix_cell{ii} > M1thresh;
%                 M2positive_master(istart:iend) =...
%                     M2pix_cell{ii} > M2thresh;
%             else
%                 M1positive_master(istart:iend) =...
%                     M1pix_cell{ii} > M1thresh_s(ii);
%                 M2positive_master(istart:iend) =...
%                     M2pix_cell{ii} > M2thresh_s(ii);
%             end
%     
%         end
   %}
   
    % Centroids
    Centroids = regionprops(Dapireg(:,:,ii),'Centroid');
    Centroids_master(istart:iend,:) = round(cell2mat({Centroids.Centroid}'));
    
    % Book keeping
    M1_cellind_master(istart:iend) = 1 : NROIs_master(ii);
    M1_slice_master(istart:iend) = ii;
end
%{
% if Automarkerthresh == 1
%     M1pix_dist = fitdist(M1pix_master,'Normal');
%     M2pix_dist = fitdist(M2pix_master,'Normal');
%
%     M1thresh_auto = M1pix_dist.mu + M1pix_dist.sigma * 1.5;
%     M2thresh_auto = M2pix_dist.mu + M2pix_dist.sigma * 0.5;
%
%     M1positive_master = M1pix_master >= M1thresh_auto;
%     M2positive_master = M2pix_master >= M1thresh_auto;
% end


%% Make centroid map

% Define brain region with a polygon
% polyroi = getpoly(max(Dapistack,[],3));

% Determine if a centroid is in the polygon or not
% inpoly_master = zeros(sum(NROIs_master,1),1);
%
% for i = 1 : sum(NROIs_master,1)
%     inpoly_master(i) = polyroi(Centroids_master(i,2), Centroids_master(i,1));
% end
%}

%% Master data matrix

% Master_data_mat = [M1_slice_master, M1_cellind_master, Centroids_master, M1pix_master,...
%     M2pix_master, inpoly_master, M1positive_master, M2positive_master];
Master_data_mat = [M1_slice_master, M1_cellind_master, Centroids_master, M1pix_master,...
    M2pix_master, M1positive_master, M2positive_master];

%{ % n_Marker1 = sum(Master_data_mat(:,7) .* Master_data_mat(:,8))
% n_Marker2 = sum(Master_data_mat(:,7) .* Master_data_mat(:,9))
% ntotal = sum(Master_data_mat(:,7))
% n_both = sum(Master_data_mat(:,7) .* Master_data_mat(:,8) .* Master_data_mat(:,9))
% chance = n_Marker1 * n_Marker2 / ntotal

%% Make plot

% scatter(mat2gray(M1pix_master), mat2gray(M2pix_master), [],... [M1positive_master*0.8, M2positive_master*0.6,zeros(sum(NROIs_master,1),1)])
%    scatter(M1pix_master, M2pix_master, [],...
%     [Master_data_mat(:,8)*0.8, Master_data_mat(:,9)*0.6,zeros(sum(NROIs_master,1),1)])
% xlabel(Marker1)
% ylabel(Marker2)

%% Slice progression

%{
F_M1pos = zeros(nslices, 1);
F_M2pos = zeros(nslices, 1);

for ii = 1 : nslices
    F_M1pos(ii) = mean(M1pix_cell{ii} > M1thresh(ii));
    F_M2pos(ii) = mean(M2pix_cell{ii} > M2thresh(ii));
end

plot(1:nslices,F_M1pos, 1:nslices,F_M2pos)
ylabel('Expression')
xlabel('slice')
legend({Marker1, Marker2})
%}

%% Map
% if Smooththresh == 1 Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh_s, M2thresh_s, NROIs_master) else Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh, M2thresh, NROIs_master) end
%}

%% Save
savefileid = fullfile(filepath,[Dapifile(1:end-3),'mat']);
save(savefileid);


%% Backup
%{ 
%%Stephen's backup

% for   ii = 1 : 3
% %     waitbar(ii/nslices)
%
%     % Image to cleanup
%     img2clean = Dapireg(:,:,ii);
%
%     % If the pixel is positive in previous three slices
%     pixel_before = max(Dapireg(:,:,max(1,(ii-3)) : max(1,(ii-1))),[],3);
%     % If the pixel is positive in the three following slices
%     pixel_after = max(Dapireg(:,:,min(ii+1,nslices) : min(ii+3,nslices)),[],3);
%     % If the pixel is positive only in current slice, mark 1, otherwise, 0 (logical).
%     tempoverlap = ((Dapireg(:,:,ii)>0) - (pixel_before>0) - (pixel_after>0))>0;
%     % Remove the pixels also positve in other slices, only keep the unique pixels in the slices.
%     tempoverlap_l = Dapireg(:,:,ii) .* tempoverlap;
%
%     % Calculate how many cells are there in current slice
%     n_areas = max(max(Dapireg(:,:,ii)));
%     % Prime parameters
%     areas_before = zeros(n_areas , 1);
%     areas_after = zeros(n_areas , 1);
%     % Examine all the cells one by one
%     for i = 1 : n_areas
%         % areas_before, how many pixels Cell i has in current slice.
%         areas_before(i) = sum(sum(Dapireg(:,:,ii) == i));
%         % area_after, how many pixels Cell i still have after subtraction (unique pixels in current slices.
%         areas_after(i) = sum(sum(tempoverlap_l == i));
%         % If current cell is not the biggest, area_after(i) == 0.
%         % ALso remove Cells with big portion left,
%         % usually because two or more cells are recognized as single cells in some but not all slices.
%         if   areas_after(i)/areas_before(i) <= 0 || areas_after(i)/areas_before(i) > 0.2
%             img2clean(img2clean==i) = 0;
%         end
%     end
%
%        figure; imshow(img2clean)
%         % figure;imshowpair(Dapistack_scaled(:,:,ii),bwperim(img2clean));
%     % Relabel and store
%     Dapireg(:,:,ii) = bwlabel(img2clean>0);
%     %figure; imshow(Dapireg(:,:,ii))
% end
    
% % If the cell is also exist in the next slice, remove it
% for ii = 2 : 3
%     img2clean = Dapireg(:,:,ii);
%     % Remove the signal exist in the next slices
%     tempoverlap = ((Dapireg(:,:,ii)>0) - (Dapireg(:,:,max(ii+1,nslices))>0))>0;
%     tempoverlap_l = Dapireg(:,:,ii) .* tempoverlap;
%     n_areas = max(max(Dapireg(:,:,ii)));
%     areas_current = zeros(n_areas , 1);
%     areas_subtracted = zeros(n_areas , 1);
%     for i = 1:n_areas
%         areas_current(i) = sum(sum(Dapireg(:,:,ii) == i));
%         areas_subtracted(i) = sum(sum(tempoverlap_l == i));
%         if areas_subtracted(i)/areas_current(i) < 0.8
%             img2clean(img2clean==i) = 0;
%         end
%     end
%     Dapireg(:,:,ii) = bwlabel(img2clean>0);
%     figure; imshow(Dapireg(:,:,ii))
% end

%     % Calculate and applyed labels to overlap
%     % Only those points with new signal comparing with the previous 5 slices will be 1.
%      tempoverlap = ((Dapireg(:,:,ii)>0) ...
%          - (max(Dapireg(:,:,max(1,(ii-overlapretro)) : (ii-1)),[],3)>0))>0;
%     tempoverlap_l = Dapireg(:,:,ii) .* tempoverlap;
%     % figure; imshow(tempoverlap);
%
%     % Calculate how many cells are there in current slice
%     n_areas = max(max(Dapireg(:,:,ii)));
%     areas_before = zeros(n_areas , 1);
%     areas_after = zeros(n_areas , 1);
%
%     % Loop through all the regions
%     for i = 1 : n_areas
%         % areas_before, how many pixels Cell i has in current slice.
%         % area_after, how many pixels Cell i still have after ...
%         % subtracting the pixels positive in the previous 5 slices.
%         areas_before(i) = sum(sum(Dapireg(:,:,ii) == i));
%         areas_after(i) = sum(sum(tempoverlap_l == i));
%
%         if areas_after(i)/areas_before(i) < (1 - overlaptreshold)
%             img2clean(img2clean==i) = 0;
%         end
%     end
%

% end

% close(hwait)
%
% for   ii = 1 : nslices
%     % Image to cleanup
%     img2clean = Dapireg(:,:,ii);
%     % Calculate how many cells are there in current slice
%     n_areas = max(max(Dapireg(:,:,ii)));
%     % Prime parameters
%     areas_current = zeros(n_areas , 1);
%     areas_overlap = zeros(n_areas , 1);
%     % Label the pixels both positive in current and next slice
%     tempoverlap = ((Dapireg(:,:,ii)>0) + (Dapireg(:,:,min(ii+1,nslices))>0))>1;
%     % Assign cell numbers of current slices to the overlap slices
%     tempoverlap_l = Dapireg(:,:,ii) .* tempoverlap;
%
%     % Examine all the cells one by one
%     for i = 1 : n_areas
%         areas_current(i) = sum(sum(Dapireg(:,:,ii) == i));
%         areas_overlap(i) = sum(sum(tempoverlap_l == i));
%         % If the overlap part is larger then 0.82 of the cell, remove this cell
%         if   areas_overlap(i)/areas_current(i) > 0.82 || areas_current(i) < 40
%             img2clean(img2clean==i) = 0;
%         end
%     end
%     % figure; imshow(img2clean);
%
%     % Relabel and store
%     Dapireg(:,:,ii) = bwlabel(img2clean>0);
%     figure; imshow(Dapireg(:,:,ii))
% end
%
%     % For testing, to examine overlaps between slices
%     %  figure;imshowpair(imadjust(Dapireg(:,:,1)),imadjust(Dapireg(:,:,2)));
%     %   figure;imshowpair(imadjust(Dapireg(:,:,12)),imadjust(Dapireg(:,:,14)));
%     % figure;imshowpair(Dapistack_scaled(:,:,12),bwperim(Dapireg(:,:,12)));
%     % figure;imshowpair(Dapistack_scaled(:,:,13),bwperim(Dapireg(:,:,13)));

%}

