%% Parameters

Smooththresh = 1; % Determine if smooth threshold or not
Automarkerthresh = 0; % If automatedly determine marker threshold

firstslice2read = 3; % First slice to read
maxnslices = 15; % Only consider first X slides.
openclose = 5; % Size to use for imopen
backgroundsize = 10; % Size to use for background subtraction
overlaptreshold = 0.5; % Fraction overlap to be considered the same cell
overlapretro = 10; % Use the past N slices for overlap considerations

Marker1 = 'LMO';
Marker2 = 'LHX1';

Commonfilapath = 'D:\Dropbox\Crickmore_research\Images\Lim1\mouse brain 3-full SCN';

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

%% Segmentate Dapi

hwait = waitbar(0, 'Segmentating');

for ii = 1 : nslices
    waitbar(ii/nslices)
    
    % Calculate background for dapi
    % dapibg = imopen(Dapistack(:,:,ii),strel('disk',35));
    
    % Logistic transform the pixel values
    % dapilog = log(1+double(mat2gray(Dapistack(:,:,ii))));
    
    % Open connections
    dapiopen = imopen(Dapistack(:,:,ii), strel('disk', openclose));
    
    % Calculate a threshold
    [thresholdvalue, areamin] = thresholdgradient(dapiopen, [200, 700]);
    
    % Apply thresholds
    dapiroi = dapiopen >= thresholdvalue;
    
    % Watershed
    watershedded = watershed( -bwdist(~dapiroi));
    dapiroi_ws = dapiroi;
    dapiroi_ws(watershedded == 0) = 0;
        
    % Apply area threshold
    dapiroi_ws_ath = areathresh( dapiroi_ws, 0.5, areamin, 3, 8);
    
    % Store the registered image
    Dapireg(:,:,ii) = bwlabel(dapiroi_ws_ath);
end

close(hwait)
%% Clean up overlaps

hwait = waitbar(0, 'Cleaning up overlap');

for ii = 2 : nslices
    waitbar(ii/nslices)
    
    % Image to cleanup
    img2clean = Dapireg(:,:,ii);

    % Calculate and applyed labels to overlap
    tempoverlap = ((Dapireg(:,:,ii)>0) ...
        - (max(Dapireg(:,:,max(1,(ii-overlapretro)) : (ii-1)),[],3)>0))>0;
    tempoverlap_l = Dapireg(:,:,ii) .* tempoverlap;

    % Calculate before and after areas
    n_areas = max(max(Dapireg(:,:,ii)));
    areas_before = zeros(n_areas , 1);
    areas_after = zeros(n_areas , 1);

    % Loop through all the regions
    for i = 1 : n_areas
        % Calculate the areas
        areas_before(i) = sum(sum(Dapireg(:,:,ii) == i));
        areas_after(i) = sum(sum(tempoverlap_l == i));

        if areas_after(i)/areas_before(i) < (1 - overlaptreshold)
            img2clean(img2clean==i) = 0;
        end
    end
    
    % Relabel and store
    Dapireg(:,:,ii) = bwlabel(img2clean>0);

end

close(hwait)

%% Prime pixel values for markers

% Prime a threshold vector
M1thresh = zeros(nslices, 1);
M2thresh = zeros(nslices, 1);

% Prime cells to contain pixel values
M1pix_cell = cell(nslices, 1);
M2pix_cell = cell(nslices, 1);

%% Calculate pixel values for markers
for ii = 1 : nslices

    % flatten backgrounds
    M1bg = imopen(Marker1stack(:,:,ii), strel('disk', backgroundsize));
    M2bg = imopen(Marker2stack(:,:,ii), strel('disk', backgroundsize));

    M1_nobg = Marker1stack(:,:,ii) - M1bg;
    M2_nobg = Marker2stack(:,:,ii) - M2bg;

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

    % Establish pixel intensity order
    [M1pix_sorted , M1order] = sort(M1pix, 1, 'descend');
    [M2pix_sorted , M2order] = sort(M2pix, 1, 'descend');
    
    % Use antomated or manual ways to determine threshold
    if Automarkerthresh == 1
        % mean + 2 std
        % M1thresh(ii) = mean(M1_nobg(:)) + 1.5 * std(M1_nobg(:));
        % M2thresh(ii) = mean(M1_nobg(:)) + 1.5 * std(M2_nobg(:));;
        
    else
        % Manual
        if M1thresh(ii) <= 0
            figurename = ['Slice ', num2str(ii), ' - Marker 1'];

            % Manually find threshold if needed
            if max(max(Dapireg(:,:,ii))) > 0
                markerthreshold(Marker1stack(:,:,ii), Dapireg(:,:,ii),...
                    M1order, figurename);
                uiwait()
                M1thresh(ii) = M1pix_sorted(lastpositive);
            else
                M1thresh(ii) = NaN;
            end


        end

        if M2thresh(ii) <= 0
            figurename = ['Slice ', num2str(ii), ' - Marker 2'];

            % Manually find threshold if needed
            if max(max(Dapireg(:,:,ii))) > 0
                markerthreshold(Marker2stack(:,:,ii), Dapireg(:,:,ii),...
                    M2order, figurename);
                uiwait()
                M2thresh(ii) = M2pix_sorted(lastpositive);
            else
                M1thresh(ii) = NaN;
            end

        end
    end

    % Load pixel values into cells
    M1pix_cell{ii} = M1pix;
    M2pix_cell{ii} = M2pix;

end


%% Smooth threshold
if Smooththresh == 1
    M1thresh_s = smooth(M1thresh);
    M2thresh_s = smooth(M2thresh);
end

%% Consolidate data


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

% Loop throug the slices
for ii = 1 : nslices
    % Determine start and end indices
    if ii == 1
        istart = 1;
    else
        istart = 1 + sum(NROIs_master(1 : (ii-1)));
    end
    
    iend = sum(NROIs_master(1:ii));
    
    % Load the values
    M1pix_master(istart:iend) = M1pix_cell{ii};
    M2pix_master(istart:iend) = M2pix_cell{ii};
    
    % Determine positivity
    if Automarkerthresh < 1
        if Smooththresh < 1
            M1positive_master(istart:iend) =...
                M1pix_cell{ii} >= M1thresh(ii);
            M2positive_master(istart:iend) =...
                M2pix_cell{ii} >= M2thresh(ii);
        else
            M1positive_master(istart:iend) =...
                M1pix_cell{ii} >= M1thresh_s(ii);
            M2positive_master(istart:iend) =...
                M2pix_cell{ii} >= M2thresh_s(ii);
        end
      
    end
    
    % Centroids
    Centroids = regionprops(Dapireg(:,:,ii),'Centroid');
    Centroids_master(istart:iend,:) = round(cell2mat({Centroids.Centroid}'));
    
    % Book keeping
    M1_cellind_master(istart:iend) = 1 : NROIs_master(ii);
    M1_slice_master(istart:iend) = ii;
end

if Automarkerthresh == 1
    M1pix_dist = fitdist(M1pix_master,'Normal');
    M2pix_dist = fitdist(M2pix_master,'Normal');

    M1thresh_auto = M1pix_dist.mu + M1pix_dist.sigma * 1.5;
    M2thresh_auto = M2pix_dist.mu + M2pix_dist.sigma * 0.5;

    M1positive_master = M1pix_master >= M1thresh_auto;
    M2positive_master = M2pix_master >= M1thresh_auto;
end


%% Make centroid map

% Define brain region with a polygon
polyroi = getpoly(max(Dapistack,[],3));

% Determine if a centroid is in the polygon or not
inpoly_master = zeros(sum(NROIs_master,1),1);

for i = 1 : sum(NROIs_master,1)
    inpoly_master(i) = polyroi(Centroids_master(i,2), Centroids_master(i,1));    
end

%% Master data matrix
Master_data_mat = [M1_slice_master, M1_cellind_master, Centroids_master, M1pix_master,...
    M2pix_master, inpoly_master, M1positive_master, M2positive_master];

%% Make plot
% scatter(mat2gray(M1pix_master), mat2gray(M2pix_master), [],...
%     [M1positive_master*0.8, M2positive_master*0.6,zeros(sum(NROIs_master,1),1)])
scatter(M1pix_master, M2pix_master, [],...
    [M1positive_master*0.8, M2positive_master*0.6,zeros(sum(NROIs_master,1),1)])
xlabel(Marker1)
ylabel(Marker2)

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

if Smooththresh == 1 
    Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh_s, M2thresh_s, NROIs_master)
else
    Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh, M2thresh, NROIs_master)
end

%% Save
savefileid = fullfile(filepath,[Dapifile(1:end-3),'mat'])
save(savefileid);