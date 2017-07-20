function imwriteallraw(filename,img,precision)

% '*uint8' 8 bit imaging, raw data, behavioral camera
% '*uint16' 16 bit imaging, raw data, VSD camera
% '*float32' 32 bit, filtered data, VSD camera

img = permute(img, [2, 1, 3]);
fid = fopen(filename,'w', 'b');
fwrite(fid, img, precision);
fclose(fid);
fprintf(1,'%s\n',filename);

end