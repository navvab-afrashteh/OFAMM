function ImgSeq = addNoiseToImgSeq(ImgSeq, NoverSig)

dim1 = size(ImgSeq,1);
dim2 = size(ImgSeq,2);
dim3 = size(ImgSeq,3);
II = ImgSeq .* ImgSeq;
sig = sum(II(:)/(dim1*dim2*dim3));
nLevel = NoverSig * sig / 100;
noise = randn(dim1, dim2, dim3) * sqrt(nLevel);
ImgSeq = ImgSeq + noise;
