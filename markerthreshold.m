function [  ] = markerthreshold( Markerstack, ROIstack, ROIrank, text )
%MARKERTHRESHOLD lets user to use a slider to determine the threshold of
%the current marker
%   [ threshold_out ] = markerthreshold( Markerstack, ROIstack, ROIrank, text )


% No text situation
if nargin < 4
    text = 'Current channel';
end

% Current rank to show
rank2show = 1;

% Make RGB
im2show = repmat(mat2gray(Markerstack), [1 1 3]);
im2show(:,:,3) = 0;
im2show(:,:,1) = (ROIstack == ROIrank(rank2show))*0.5;

figure('Position', [50, 50, 1400, 700])
h1 = imshow(im2show, []);

set(gcf,'Name', text);

% Add slider
hsl = uicontrol('Style', 'Slider','Position', [50 20 320 20],...
    'Min',1, 'Max', length(ROIrank),'Value',1,...
    'SliderStep',[1/length(ROIrank) 10/length(ROIrank)],...
    'Callback', @changeroi);

% Add a button
hbt =  uicontrol('Style', 'pushbutton', 'String', 'Confirm',...
    'Position', [400, 20, 80, 20], 'Callback', @confirmcurrent);

    function changeroi(source,~)
        im2show(:,:,1) = (ROIstack == ROIrank(round(source.Value)))*0.5;
        h1.CData = im2show;
    end

    function confirmcurrent(~,~)
        currentchoice = round(hsl.Value);
        assignin('base','lastpositive', currentchoice)
        close
    end


end

