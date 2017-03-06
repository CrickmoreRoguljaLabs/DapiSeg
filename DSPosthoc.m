
Commonfilapath = 'D:\Dropbox\Crickmore_research\Images\Lim1\mouse brain 3-full SCN';

[matfiles, filepath2] = uigetfile(fullfile(Commonfilapath, '*.mat'), 'MultiSelect', 'on');

nfiles = length(matfiles);

Master_data_cell = cell(nfiles,1);
Master_data_mat_posthoc = [];


for iii = 1 : nfiles
    disp(iii)
    filehandle = matfiles{iii};
    load(fullfile(filepath2, filehandle));
       
    
    %Markermap( Dapireg, M1pix_cell, M2pix_cell, M1thresh_s, M2thresh_s, NROIs_master)
    %print(fullfile(filepath2, [filehandle(1:3),'_map.png']),'-dpng')
    
    if iii == 1
        Master_data_mat_posthoc = Master_data_mat;
        counter = size(Master_data_mat_posthoc, 1);
    else
        n2add = size(Master_data_mat,1);
        Master_data_mat_posthoc((counter+1) : (counter+n2add), :)  = Master_data_mat;
        counter = counter + n2add;
    end
        
    Master_data_cell{iii} = Master_data_mat;
end

keep Commonfilapath nfiles Master_data_cell Master_data_mat_posthoc matfiles filepath2
%%
M1fraction = zeros(nfiles,1);
M2fraction = zeros(nfiles,1);

for iii = 1 : nfiles
    currentstat = mean(Master_data_cell{iii});
    
    M1fraction(iii) = currentstat(1,5);
    M2fraction(iii) = currentstat(1,6);
end
%%
plot(1:nfiles,M1fraction, 1:nfiles,M2fraction)
ylabel('Marker/Dapi')
xlabel('slice')
legend({'LMO', 'Lim1'})
%%
scatter(Master_data_mat_posthoc(:,3), Master_data_mat_posthoc(:,4), [],...
    [Master_data_mat_posthoc(:,5)*0.8, Master_data_mat_posthoc(:,6)*0.6,...
    zeros(size(Master_data_mat_posthoc,1),1)])

xlabel('LMO level')
ylabel('Lim1 level')