function imwriteallraw(filename,img,precision)

    % '*uint8' 8 bit imaging, raw data, behavioral camera
    % '*uint16' 16 bit imaging, raw data, VSD camera
    % '*float32' 32 bit, filtered data, VSD camera
    
   if size(img,3) > 1
       for i = 1: size(img,3)
            img1(:,:,i) = fliplr(rot90(squeeze(img(:,:,i)),-1));         
       end
   else
       img1 = fliplr(rot90(img,-1));
   end
 
   fid = fopen(filename,'w', 'b');
   fwrite(fid, img1, precision);
   fclose(fid);
   fprintf(1,'%s\n',filename);

end