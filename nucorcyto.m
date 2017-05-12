%% Load the variables that are used in this program

[matfile, matpath] = uigetfile('D:\Dropbox\Brain 3 SCN mat files\*.mat', 'Select a workspace file from DapiSeg');
load(fullfile(matpath, matfile), 'Master_data_mat', 'Dapistack', 'Marker1stack', 'Marker2stack', 'Marker1', 'Marker2')


%% Make stack of all M2 positive cells
% How many pixels to extend from each centroid
expansionsize = 6;

% Calculate the real dimension of the stacks
expansionsize_real = expansionsize * 2 + 1;

% Find the indices of Marker 2 positive cells
M2_pos_ind = find(Master_data_mat(:,8)>0);

% The limits of pixel indices
maxY = size(Marker2stack,2);
MaxX = size(Marker2stack,1);

% Get the number of Marker 2 cells and prime the matrices containing the
% M2 and Dapi images
n_M2 = length(M2_pos_ind);
M2_pos_stack = zeros(expansionsize_real, expansionsize_real, n_M2);
M2_pos_stack_dapi = zeros(expansionsize_real, expansionsize_real, n_M2);
M2_pos_stack_M1 = zeros(expansionsize_real, expansionsize_real, n_M2);

for ii = 1 : n_M2
    
    % Locate the index of the cell
    matind = M2_pos_ind(ii);
    
    % Find the X, Y, Z location ranges of the cell
    sectionind = Master_data_mat(matind,1);
    YindU = max(1, Master_data_mat(matind,3) - expansionsize);
    YindD = min(maxY, Master_data_mat(matind,3) + expansionsize);
    XindL = max(1, Master_data_mat(matind,4) - expansionsize);
    XindR = min(MaxX, Master_data_mat(matind,4) + expansionsize);

    % Calculate the dimension and the indices in the new stack (centered as
    % much as possible(
    Ydim = YindD - YindU + 1;
    Xdim = XindR - XindL + 1;
    Xstart = floor((expansionsize_real - Xdim)/2) +  1;
    Ystart = floor((expansionsize_real - Ydim)/2) +  1;
    Xend = Xstart + Xdim - 1;
    Yend = Ystart + Ydim - 1;
    
    % Load the cell unto the new stack
    M2_pos_stack(Xstart:Xend,Ystart:Yend,ii) = Marker2stack(XindL:XindR,...
        YindU:YindD, sectionind);
    M2_pos_stack_dapi(Xstart:Xend,Ystart:Yend,ii) = Dapistack(XindL:XindR,...
        YindU:YindD, sectionind);
    M2_pos_stack_M1(Xstart:Xend,Ystart:Yend,ii) = Marker1stack(XindL:XindR,...
        YindU:YindD, sectionind);
end


% Play the stack if needed
% [M2, Dapi]
% stack2show_M2 = [M2_pos_stack, M2_pos_stack_dapi];
% [M1, M2, Dapi]
stack2show_M2 = [M2_pos_stack_M1, M2_pos_stack, M2_pos_stack_dapi];
implay(stack2show_M2/2000)

%% Make stack of all M1 positive cells
% How many pixels to extend from each centroid
expansionsize = 6;

% Calculate the real dimension of the stacks
expansionsize_real = expansionsize * 2 + 1;

% Find the indices of Marker 1 positive cells
M1_pos_ind = find(Master_data_mat(:,7)>0);

% The limits of pixel indices
maxY = size(Marker1stack,2);
MaxX = size(Marker1stack,1);

% Get the number of Marker 1 cells and prime the matrices containing the
% M2 and Dapi images
n_M1 = length(M1_pos_ind);
M1_pos_stack = zeros(expansionsize_real, expansionsize_real, n_M1);
M1_pos_stack_dapi = zeros(expansionsize_real, expansionsize_real, n_M1);
% M2_pos_stack_M1 = zeros(expansionsize_real, expansionsize_real, n_M2);

for ii = 1 : n_M1
    
    % Locate the index of the cell
    matind = M1_pos_ind(ii);
    
    % Find the X, Y, Z location ranges of the cell
    sectionind = Master_data_mat(matind,1);
    YindU = max(1, Master_data_mat(matind,3) - expansionsize);
    YindD = min(maxY, Master_data_mat(matind,3) + expansionsize);
    XindL = max(1, Master_data_mat(matind,4) - expansionsize);
    XindR = min(MaxX, Master_data_mat(matind,4) + expansionsize);

    % Calculate the dimension and the indices in the new stack (centered as
    % much as possible(
    Ydim = YindD - YindU + 1;
    Xdim = XindR - XindL + 1;
    Xstart = floor((expansionsize_real - Xdim)/2) +  1;
    Ystart = floor((expansionsize_real - Ydim)/2) +  1;
    Xend = Xstart + Xdim - 1;
    Yend = Ystart + Ydim - 1;
    
    % Load the cell unto the new stack
    M1_pos_stack(Xstart:Xend,Ystart:Yend,ii) = Marker1stack(XindL:XindR,...
        YindU:YindD, sectionind);
    M1_pos_stack_dapi(Xstart:Xend,Ystart:Yend,ii) = Dapistack(XindL:XindR,...
        YindU:YindD, sectionind);

end

% stack2show = [M2_pos_stack_M1, M2_pos_stack, M2_pos_stack_dapi];
% Play the stack if needed
stack2show_M1 = [M1_pos_stack, M1_pos_stack_dapi];
implay(stack2show_M1/3000)


%% Make co-centric curves for M2
% Prime the matricies containing 
Vecmean_M2 = zeros(n_M2, expansionsize_real);
Vecmean_M2_dapi = zeros(n_M2, expansionsize_real);
Vecmean_M2_M1 = zeros(n_M2, expansionsize_real);


for ii = 1 : n_M2
    % 4 vectors, the average of which describes the concavity of Marker1
    % fluorescence
    
    % Isolate the current iamge
    M2_pos_current = M2_pos_stack(:,:,ii);
    
    % Horizontal pixels
    Vec1 = M2_pos_current(expansionsize + 1, :);
    
    % Vertical pixels
    Vec2 = M2_pos_current(:, expansionsize + 1);
    
    % One set of diagonal pixels
    Vec3 = M2_pos_current(1 : expansionsize_real+1 : expansionsize_real^2);
    
    % Another set of diagonal pixels
    Vec4 = M2_pos_current(expansionsize_real : expansionsize_real-1 : expansionsize_real^2-1);

    % Average the four dimensions and Load the values to matrix
    Vecmean_M2(ii,:) = mat2gray(mean([Vec1; Vec2'; Vec3; Vec4]));
    
    % Repeat for Dapi
    M2_pos_current_dapi = M2_pos_stack_dapi(:,:,ii);
    Vec1_dapi = M2_pos_current_dapi(expansionsize + 1, :);
    Vec2_dapi = M2_pos_current_dapi(:, expansionsize + 1);
    Vec3_dapi = M2_pos_current_dapi(1 : expansionsize_real+1 : expansionsize_real^2);
    Vec4_dapi = M2_pos_current_dapi(expansionsize_real : expansionsize_real-1 : expansionsize_real^2-1);
    
    % Load Dapi values
    Vecmean_M2_dapi(ii,:) = mat2gray(mean([Vec1_dapi; Vec2_dapi'; Vec3_dapi; Vec4_dapi]));
    
    % Repeat for Marker 1
    M2_pos_current_M1 = M2_pos_stack_M1(:,:,ii);
    Vec1_M1 = M2_pos_current_M1(expansionsize + 1, :);
    Vec2_M1 = M2_pos_current_M1(:, expansionsize + 1);
    Vec3_M1 = M2_pos_current_M1(1 : expansionsize_real+1 : expansionsize_real^2);
    Vec4_M1 = M2_pos_current_M1(expansionsize_real : expansionsize_real-1 : expansionsize_real^2-1);
    
    % Load Dapi values
    Vecmean_M2_M1(ii,:) = mat2gray(mean([Vec1_M1; Vec2_M1'; Vec3_M1; Vec4_M1]));
end


%% Make co-centric curves for M1
% Prime the matricies containing 
Vecmean_M1 = zeros(n_M1, expansionsize_real);
Vecmean_M1_dapi = zeros(n_M1, expansionsize_real);


for ii = 1 : n_M1
    % 4 vectors, the average of which describes the concavity of Marker1
    % fluorescence

    % Isolate the current iamge
    M1_pos_current = M1_pos_stack(:,:,ii);
    
    % Horizontal pixels
    Vec1 = M1_pos_current(expansionsize + 1, :);
    
    % Vertical pixels
    Vec2 = M1_pos_current(:, expansionsize + 1);
    
    % One set of diagonal pixels
    Vec3 = M1_pos_current(1 : expansionsize_real+1 : expansionsize_real^2);
    
    % Another set of diagonal pixels
    Vec4 = M1_pos_current(expansionsize_real : expansionsize_real-1 : expansionsize_real^2-1);

    % Average the four dimensions and Load the values to matrix
    Vecmean_M1(ii,:) = mat2gray(mean([Vec1; Vec2'; Vec3; Vec4]));
    
    % Repeat for Dapi
    M1_pos_current_dapi = M1_pos_stack_dapi(:,:,ii);
    Vec1_dapi = M1_pos_current_dapi(expansionsize + 1, :);
    Vec2_dapi = M1_pos_current_dapi(:, expansionsize + 1);
    Vec3_dapi = M1_pos_current_dapi(1 : expansionsize_real+1 : expansionsize_real^2);
    Vec4_dapi = M1_pos_current_dapi(expansionsize_real : expansionsize_real-1 : expansionsize_real^2-1);

    % Load Dapi values
    Vecmean_M1_dapi(ii,:) = mat2gray(mean([Vec1_dapi; Vec2_dapi'; Vec3_dapi; Vec4_dapi]));
end


%% Page views for Marker 2 Co-centric curves
%{
% Current page number, can change
pagenum = 1;

% Define the pages
panelsperpage = 25;

% Initiate the figure
figure('Position',[20 50 1500 700])
for i = 1 : panelsperpage
    
    % Subplot
    subplot(5,5,i)
    
    % Calculate the panel ID
    panelind = 25 * (pagenum-1) + i;
    
    if panelind <= n_M2
        % Make the plot
        plot(1:expansionsize_real, Vecmean_M2(panelind,:),1:expansionsize_real, Vecmean_M2_dapi(panelind,:))

        % Lobel with the Marker 2 ID, different from the global Cell ID
        title(['Marker2 ID: ', num2str(panelind)])
        
        % Legend
        if i == panelsperpage
            legend({Marker2;'Dapi'})
        end
    end
    
end
%}
%% Page views for Marker 1 Co-centric curves
%{
% Current page number, can change
pagenum = 1;

% Define the pages
panelsperpage = 25;

% Initiate the figure
figure('Position',[20 50 1500 700])
for i = 1 : panelsperpage
    
    % Subplot
    subplot(5,5,i)
    
    % Calculate the panel ID
    panelind = 25 * (pagenum-1) + i;
    
    if panelind <= n_M1
        % Make the plot
        plot(1:expansionsize_real, Vecmean_M1(panelind,:),1:expansionsize_real, Vecmean_M1_dapi(panelind,:))

        % Lobel with the Marker 1 ID, different from the global Cell ID
        title(['Marker1 ID: ', num2str(panelind)])

        % Legend
        if i == panelsperpage
            legend({Marker1;'Dapi'})
        end
    end
    
end
%}


%% Final calculations
% Relate indices of the double positive cells
Co_pos_ind_M2 = Master_data_mat(M2_pos_ind,7); % In M2 index
Co_pos_ind_M1 = Master_data_mat(M1_pos_ind,8); % In M1 index


% Prime the vectors of the cytoplasmic cells
cytovec_M1 = zeros(n_M1,1);
cytovec_M2 = zeros(n_M2,1);

% Enter the indices here
% This part is changed by the user
cytoind_M1 = [];
cytoind_M2 = [30 36 44 45 48 49 58 60 129 137 154 155 165 174 175 191 201 213 221 304];





 
%%
% Make a figure of all Marker 2 cytoplasmic cells
figure('Position',[20 50 1500 700])
for i = 1 : length(cytoind_M2)
    % Subplot
    if i <= 25
        subplot(5,5,i)
    end
    
    % Take care of the panel ID and the vector
    panelind = cytoind_M2(i);
    cytovec_M2(panelind) = 1;
    
    if panelind <= n_M2 && i <= 25
        % Plot the curves
        if Co_pos_ind_M2(i) == 1
            % If Marker1 is positive, plot Marker1 as well
            plot(1:expansionsize_real, Vecmean_M2(panelind,:),...
                1:expansionsize_real, Vecmean_M2_dapi(panelind,:),...
                1:expansionsize_real, Vecmean_M2_M1(panelind,:))
            legend({Marker2;'Dapi';Marker1})
        else
            % If Marker1 is not positive, don't plot Marker 1
            plot(1:expansionsize_real, Vecmean_M2(panelind,:),...
                1:expansionsize_real, Vecmean_M2_dapi(panelind,:))
        end
        set(gca,'XTick',[]);
        
        % Use Marker 2 ID
        title(['Marker 2 ID: ', num2str(panelind)])

        % Legend
        if i == length(cytoind_M2)
            legend({Marker2;'Dapi'})
        end
    end
    
end

% Make a figure of all Marker 1 cytoplasmic cells
figure('Position',[20 50 1500 700])
for i = 1 : length(cytoind_M1)
    % Subplot
    subplot(5,5,i)
    
    % Take care of the panel ID and the vector
    panelind = cytoind_M1(i);
    cytovec_M1(panelind) = 1;
    
    if panelind <= n_M1
        % Plot the curves
        plot(1:expansionsize_real, Vecmean_M1(panelind,:),1:expansionsize_real, Vecmean_M1_dapi(panelind,:))
        set(gca,'XTick',[]);

        % Use Marker1 ID
        title(['Marker 1 ID: ', num2str(panelind)])

        % Legend
        if i == length(cytoind_M1)
            legend({Marker1;'Dapi'})
        end
    end
    
end

%% Calculations
% Display the calculations
disp('Total number of cells:')
sum(Master_data_mat(:,5)>0)

disp('Total number of LMO3 cells:')
sum(Master_data_mat(:,7)>0)

disp('Total number of LHX1 cells:')
sum(Master_data_mat(:,8)>0)

disp('Total number of double-positive cells:')
sum(Co_pos_ind_M2>0)

disp('Total number of cytoplasmic-LHX1 cells:')
length(cytoind_M2)

disp('Total number of double-positive and cytoplasmic-LHX1 cells:')
sum(cytovec_M2.*Co_pos_ind_M2)
% find((cytovec_M2.*Co_pos_ind_M2)>0)

disp('Total number of cytoplasmic-LM3 cells:')
length(cytoind_M1)

disp('Total number of double-positive and cytoplasmic-LMO3 cells:')
sum(cytovec_M1.*Co_pos_ind_M1)
% find((cytovec_M1.*Co_pos_ind_M1)>0)

%%
%{
TempM2 = [Master_data_mat(M2_pos_ind,:),cytovec_M2];
TempM1 = [Master_data_mat(M1_pos_ind,:),cytovec_M1];


A=TempM1;
size(A,1);
B=A(A(:,8)==1,:);
C=A(A(:,8)==0,:);

[~,Brank] = sort(B(:,6),'descend');
Bsorted = B(Brank,:);
quartilesize = round(size(B,1)/4);
B1 = mean(Bsorted(1:quartilesize,9));
B2 = mean(Bsorted((quartilesize+1) : (quartilesize * 2),9));
B3 = mean(Bsorted((2*quartilesize+1) : (quartilesize * 3),9));
B4 = mean(Bsorted((3*quartilesize+1) : end,9));
C1 = mean(C(:,9));

[B1; B2; B3; B4; C1]
%}