
Commonfilapath = 'D:\Dropbox\Crickmore_research\Images\Lim1\mouse brain 3-full SCN';

[matfiles, filepath2] = uigetfile(fullfile(Commonfilapath, '*.mat'), 'MultiSelect', 'on');

nfiles = length(matfiles);

Master_data_cell = cell(nfiles,1);
Master_data_mat_posthoc = [];


for iii = 1 : nfiles
    disp(iii)
    filehandle = matfiles{iii};
    load(fullfile(filepath2, filehandle));
       
    inroi_ind = Master_data_mat(:,7) == 1;
    
    %Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh_s, M2thresh_s, NROIs_master)
    %print(fullfile(filepath2, [filehandle(1:3),'_map.png']),'-dpng')
    
    if iii == 1
        Master_data_mat_posthoc = Master_data_mat(inroi_ind,:);
        counter = size(Master_data_mat_posthoc, 1);
    else
        n2add = sum(inroi_ind);
        Master_data_mat_posthoc((counter+1) : (counter+n2add), :)  = ...
            Master_data_mat(inroi_ind,:);
        counter = counter + n2add;
    end
        
    Master_data_cell{iii} = Master_data_mat(inroi_ind,:);
end

keep Commonfilapath nfiles Master_data_cell Master_data_mat_posthoc matfiles filepath2
%%
M1fraction = zeros(nfiles,1);
M2fraction = zeros(nfiles,1);

for iii = 1 : nfiles
    currentstat = mean(Master_data_cell{iii});
    
    M1fraction(iii) = currentstat(1,8);
    M2fraction(iii) = currentstat(1,9);
end
%% Plot slice progression
plot(1:nfiles,M1fraction, 1:nfiles,M2fraction)
ylabel('Marker/Dapi')
xlabel('slice')
legend({'LMO', 'Lim1'})

% Plot scattering
figure
scatter(mat2gray(Master_data_mat_posthoc(:,5)),...
    mat2gray(Master_data_mat_posthoc(:,6)), [],...
    [Master_data_mat_posthoc(:,8)*0.8, Master_data_mat_posthoc(:,9)*0.6,...
    zeros(size(Master_data_mat_posthoc,1),1)])

xlabel('LMO level')
ylabel('Lim1 level')

% Plot distribution
% LMO
figure
M1norm = mat2gray(Master_data_mat_posthoc(:,5));

hist(M1norm,500)
M1mean = mean(M1norm);
M1std = std(M1norm);
hold on
plot([M1mean+M1std,M1mean+M1std],[0 90],'r')
plot([M1mean+2*M1std,M1mean+2*M1std],[0 90],'r')
hold off
ylabel('Count')
xlabel('Normalized Intensity')
title('LMO3')

% LHX1
figure
M2norm = mat2gray(Master_data_mat_posthoc(:,6));

hist(M2norm,500)
M2mean = mean(M2norm);
M2std = std(M2norm);
hold on
plot([M2mean+M2std,M2mean+M2std],[0 90],'r')
plot([M2mean+2*M2std,M2mean+2*M2std],[0 90],'r')
hold off
ylabel('Count')
xlabel('Normalized Intensity')
title('LHX1')

%%
for iii = 1 :3: 9
slice2map = iii;
mapmat = Master_data_cell{slice2map};
figure
title(num2str(iii))
scatter(mapmat(:,3),mapmat(:,4),50,1-[mapmat(:,9), mapmat(:,8),ones(size(mapmat,1),1)*0.2],'filled')
end