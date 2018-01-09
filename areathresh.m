function [ outputmat ] = areathresh( inputmat, intensitythreshold,areathreshold, dim , fouroreight)
%AREATHRESH applies a pixel-value threshold and an area-size threshold
%   [ outputmat ] = areathresh( inputmat, intensitythreshold,areathreshold, dim , fouroreight)
if nargin<5
    fouroreight=8;
end

outputmat=inputmat;

if dim==1
    for i=1:size(inputmat,dim)
        tmp_area=squeeze(abs(inputmat(i,:,:)))>intensitythreshold;
        tmp_area_labeled=uint16(bwlabel(tmp_area,fouroreight));
        tmp_area_struct=regionprops(tmp_area_labeled,'Area');
        temp_area_area=cell2mat({tmp_area_struct.Area});
        areas2eliminate=find(temp_area_area<areathreshold);
        for j=1:length(areas2eliminate)
            tmp_area(tmp_area_labeled==areas2eliminate(j))=0;
        end
        tmp_area=imclose(tmp_area,strel('Diamond',1));
        outputmat(i,:,:)=tmp_area;
    end
elseif dim==3
    for i=1:size(inputmat,dim)
        tmp_area=squeeze(abs(inputmat(:,:,i)))>intensitythreshold;
        tmp_area_labeled=uint16(bwlabel(tmp_area,fouroreight));
        tmp_area_struct=regionprops(tmp_area_labeled,'Area');
        temp_area_area=cell2mat({tmp_area_struct.Area});
        areas2eliminate=find(temp_area_area<areathreshold);
        for j=1:length(areas2eliminate)
            tmp_area(tmp_area_labeled==areas2eliminate(j))=0;
        end
        outputmat(:,:,i)=tmp_area;
    end
elseif dim==0
    tmp_area=abs(inputmat(:,:))>intensitythreshold;
    tmp_area_labeled=uint16(bwlabel(tmp_area,fouroreight));
    tmp_area_struct=regionprops(tmp_area_labeled,'Area');
    temp_area_area=cell2mat({tmp_area_struct.Area});
    areas2eliminate=find(temp_area_area<areathreshold);
    for j=1:length(areas2eliminate)
            tmp_area(tmp_area_labeled==areas2eliminate(j))=0;
    end
    outputmat(:,:)=tmp_area;
end

end

