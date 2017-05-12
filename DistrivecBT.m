function [ AmpOut, XOut, YOut ] = DistrivecBT( Xall,Yall, nUnit, nIt )
%DistrivecBT bootstrap relative distribution vectors
%   DistrivecBT( Xall,Yall, nUnit, nIt )

% Default 10000 iteractions
if nargin < 4
    nIt = 10000;
end

% Initiate the outputs
XOut = zeros(nIt, 1);
YOut = zeros(nIt, 1);

% Generate a random matrix
randmat = randi(length(Xall), nUnit, nIt);

% Relativize the coordinates
Xall2 = Xall - mean(Xall);
Yall2 = Yall - mean(Yall);

for i = 1 : nIt
    % Bootstrap the coordinates
    XOut(i) = mean(Xall2(randmat(:,i)));
    YOut(i) = mean(Yall2(randmat(:,i)));
end

% Calculate the amplitudes
AmpOut = sqrt(XOut .^ 2 + YOut .^ 2); 

end

