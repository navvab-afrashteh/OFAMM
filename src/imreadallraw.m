function ImgSeq = imreadallraw(filename,x,y,nFrames,precision)

    % '*uint8' 8 bit imaging, raw data, behavioral camera
    % '*uint16' 16 bit imaging, raw data, VSD camera
    % '*float32' 32 bit, filtered data, VSD camera
    
    fid0 = fopen(filename, 'r', 'b');
    ImgSeq = fread(fid0,[x*y nFrames],precision);
    fclose(fid0);
  
    ImgSeq = reshape(ImgSeq,x,y,nFrames);
    for ii=1:nFrames
        ImgSeq(:,:,ii) = flipud(rot90(ImgSeq(:,:,ii)));
    end
end
