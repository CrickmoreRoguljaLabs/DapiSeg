function [ outputim ] = getpoly( inputim , text)
%getpoly returns a polygon region of interest selected by the user
%   output image = getpoly (input image in single channel , title of the image)

if nargin < 2
    % Default text
    text = 'select a polygon';
end

% Show the image
figure( 101 ), imshow( inputim , [] );
set(101,'Position',[50 50 800 1600])


% Set the image title
set( 101 , 'Name' , text );

% Place polygon interactively.
hpoly = impoly(gca);

% Wait until a doubleclick confirmation
wait( hpoly );

% Output the polygon 
outputim = createMask( hpoly );

% Close the image
close(101)
end

