function [fx, fy, ft] = computeDerivatives(im1, im2)

if size(im2,1)==0
    im2=zeros(size(im1));
end

% Horn-Schunck original method
% fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
% fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
% ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');

% derivatives as in Barron
I = im1+im2;
fx= conv2(I/2,(1/12)*[-1 8 0 -8 1],'same');
fy= conv2(I/2,(1/12)*[-1 8 0 -8 1]','same');
ft = conv2(im1, 0.25*ones(2),'same') + conv2(im2, -0.25*ones(2),'same');
fx=-fx;fy=-fy;

% An alternative way to compute the spatiotemporal derivatives is to use simple finite difference masks.
% fx = conv2(im1,[1 -1], 'same');
% fy = conv2(im1,[1; -1]','same');
% ft= im2-im1;

% Sobel filter
% fx = conv2(im1,[-1 -2 -1; 0 0 0; 1 2 1], 'same');
% fy = conv2(im1,[-1 -2 -1; 0 0 0; 1 2 1]', 'same');
% ft= -im2+im1;