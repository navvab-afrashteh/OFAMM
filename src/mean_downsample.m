function [gridX, gridY, downU] = mean_downsample(U,Dsize)

if (ndims(U)==2)
    [dimX,dimY] = size(U);
    dimZ = 1;
    downU = zeros(dimX,dimY);
elseif (ndims(U)==3)
    [dimX,dimY,dimZ] = size(U);
    downU = zeros(dimX,dimY,dimZ);
else
    error('ERROR:DIMENSION','Function currently requires a 2 or 3 dimensional input array, e.g., an image or an image stack');
end
    
% Build 2D mean filter
Kernel = ones(Dsize)/(Dsize*Dsize);
% Downsample initial reference point (Center,Center)
Center = ceil(Dsize/2);
% Downsample grid
gridX  = Center:Dsize:dimX;
gridY  = Center:Dsize:dimY;
% Apply 2D mean filter
for i=1:dimZ
    downU(:,:,i) = conv2(U(:,:,i),Kernel,'same');
end
% Downsample filtered imput array
downU = downU(gridX,gridY,:);

end