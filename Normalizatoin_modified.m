% Copy Master_data_set
temp = Master_data_mat;
temp_normalized = temp;

% Extract LMO and LHX1 data, organize the data by slices (2 ~ nslices-1).
% for ii = 2 : nslices-1
%     % Get the index of the lines from the same slice
%     ind = temp(:,1) == ii;
%     % Extract the intensity data from the same slice
%     lmo_intensity = temp(ind, 5);
%     lhx_intensity = temp(ind, 6);
%     % Get the maximun and minimum intensities in the same slice
%     lmo_max = max(lmo_intensity);
%     lmo_min = min(lmo_intensity(lmo_intensity>40));
%     lhx_max = max(lhx_intensity);
%     lhx_min = min(lhx_intensity(lhx_intensity>40));
%     % Normalize the intensities
%     temp_normalized(ind,5) = (lmo_intensity - lmo_min)/(lmo_max);
%     temp_normalized(ind,6) = (lhx_intensity - lhx_min)/(lhx_max);
% end

lmo_intensity = temp(:,5);
lhx_intensity = temp(:,6);
lmo_max = max(lmo_intensity);
    lmo_min = min(lmo_intensity(lmo_intensity>40));
    lhx_max = max(lhx_intensity);
    lhx_min = min(lmo_intensity(lhx_intensity>40));
 temp_normalized(:,5) = (lmo_intensity - lmo_min)/(lmo_max);
 temp_normalized(:,6) = (lhx_intensity - lhx_min)/(lhx_max);

% Delete the data of the first and the last slices
% Decide the row number of the new set of data, start from slice 2
slice1 = find(temp(:,1) ==1, 1, 'last' ) + 1;
slice_last = find(temp(:,1) == nslices,1) - 1;
A_Normalized = temp_normalized(slice1:slice_last,:);

% Get the data of LMO-LHX1 double positive cells
ind_lmo_lhx1 = A_Normalized(:,7) == 1 & A_Normalized(:,8) == 1;
A_LMO_LHX1 = A_Normalized(ind_lmo_lhx1,:);

% Get the data of LMO-only cells
ind_lmo_lhx1 = A_Normalized(:,7) == 1 & A_Normalized(:,8) == 0;
A_LMO = A_Normalized(ind_lmo_lhx1,:);

% Get the data of LHX1-only cells
ind_lmo_lhx1 = A_Normalized(:,7) == 0 & A_Normalized(:,8) == 1;
A_LHX1 = A_Normalized(ind_lmo_lhx1,:);

% Get the data of double negative cells
ind_lmo_lhx1 = A_Normalized(:,7) == 0 & A_Normalized(:,8) == 0;
A_null = A_Normalized(ind_lmo_lhx1,:);




