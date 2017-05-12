%% Load the variables that are used in this program

% Individual plot or not
indiplot = 1;

% Obtain the matrix
[matfile, matpath] = uigetfile('D:\Dropbox\Brain 3 SCN mat files\*.mat', 'Select a workspace file from DapiSeg', 'MultiSelect', 'on');
if iscell(matfile)
    nfiles = length(matfile);
else
    nfiles = 1;
end

% Initiate a matrix containing all the vectors
M1vecmat = zeros(nfiles,2);
M2vecmat = zeros(nfiles,2);

for i = 1 : nfiles
    if iscell(matfile)
        load(fullfile(matpath, matfile{i}), 'Master_data_mat', 'Dapistack', 'Marker1', 'Marker2')
    else
        load(fullfile(matpath, matfile), 'Master_data_mat', 'Dapistack', 'Marker1', 'Marker2')
    end

    % Obtain the indices
    Indall = Master_data_mat(:,5) > 0;
    Master_data_mat2 = Master_data_mat(Indall,:);
    nCells = sum(Indall);

    % Get all X and Y values
    Xall = Master_data_mat2(:, 3);
    Yall = -Master_data_mat2(:, 4);

    % X centroid
    XCen = mean(Xall);
    YCen = mean(Yall);

    % Get the relative X Y values
    Xrel = Xall - XCen;
    Yrel = Yall - YCen;

    % Get M1 and M2 indices
    IndM1 = Master_data_mat2(:,7) == 1;
    IndM2 = Master_data_mat2(:,8) == 1;

    % Mean vectors
    M1vec  = [mean(Xrel(IndM1)), mean(Yrel(IndM1))];
    M2vec  = [mean(Xrel(IndM2)), mean(Yrel(IndM2))];
    
    if indiplot == 1
        % Make plot for each slice
        figure(i)
        scatter(Xall,Yall,[],[Master_data_mat2(:,7:8)*0.9, 0.8 * ones(nCells, 1)]);
        hold on
        quiver(XCen,YCen,M1vec(1),M1vec(2),'LineWidth',2);
        quiver(XCen,YCen,M2vec(1),M2vec(2),'LineWidth',2);
        hold off
        legend({'Cells';Marker1; Marker2})

        % figure
        % imshow(max(Dapistack,[],3),[])
    end
    % Load the values in the output matrix
    M1vecmat(i,:) = M1vec;
    M2vecmat(i,:) = M2vec;
    
end

%% Conversion factors (convert from pixel to um)
% There are probably outliers and you should write the values in;
scaling_factor = ones(nfiles,2) / 1.6103;
M1vecmat_scaled = M1vecmat .* scaling_factor;
M2vecmat_scaled = M2vecmat .* scaling_factor;

%% Make Vector plots
figure
hold on
quiver(zeros(nfiles,1),zeros(nfiles,1),M1vecmat_scaled(:,1),M1vecmat_scaled(:,2),'b')
quiver(zeros(nfiles,1),zeros(nfiles,1),M2vecmat_scaled(:,1),M2vecmat_scaled(:,2),'r')

quiver(0,0,mean(M1vecmat_scaled(:,1)),mean(M1vecmat_scaled(:,2)),'b','LineWidth',4)
quiver(0,0,mean(M2vecmat_scaled(:,1)),mean(M2vecmat_scaled(:,2)),'r','LineWidth',4)
xlabel('Lateral -> Medial (um)')
ylabel('Ventral -> Dorsal (um)')
hold off