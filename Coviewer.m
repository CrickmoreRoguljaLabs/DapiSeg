plane2view = 8;

planeind = Master_data_mat(:,1) == plane2view;
datamat_tmp = Master_data_mat(planeind,:);

colabid = datamat_tmp(:,5).* datamat_tmp(:,6);
datamat_tmp2 = datamat_tmp(colabid==1,:);

rgb2show = repmat(mat2gray(Marker1stack(:,:,plane2view)),[1 1 3]);
rgb2show(:,:,2) = mat2gray(Marker2stack(:,:,plane2view));
rgb2show(:,:,3) = 0;

if ~isempty(datamat_tmp2)
    for i = 1 : size(datamat_tmp2,1)
        rgb2show(:,:,3) = rgb2show(:,:,3) + (Dapireg(:,:,plane2view) == datamat_tmp2(i,2));
    end
end

rgb2show(:,:,3) = rgb2show(:,:,3) * 0.3;

imshow(rgb2show)
figure
imshow(rgb2show(:,:,1))
figure
imshow(rgb2show(:,:,2))
imtool(rgb2show(:,:,3)/0.3.*Dapireg(:,:,plane2view))

%%

