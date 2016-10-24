function [ motionMat ] = BatchField( imgStack )
%BATCHFIELD Summary of this function goes here
%   Detailed explanation goes here
fname = imgStack;
info = imfinfo(fname);
num_images = numel(info);
for k = 1:(num_images-1)
    A = imread(fname, k);
    B = imread(fname,k+1);
    [vx(:,:,k),vy(:,:,k),warpI2]=Coarse2FineTwoFrames(A,B);
end
motionMat=vx+i*vy;
end

