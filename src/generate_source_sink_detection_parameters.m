function generate_source_sink_detection_parameters(SSParamsPathName)

Nmin = 5; % minimum number of points as the contour size
Lmin_source = 0.05; % minimum source level
Lmax_sink = -0.05; % maximum sink level
Nnested_source = 2; % minimum number of nested sources to verify the most interior source
Nnested_sink = 2; % minimum number of nested sources to verify the most interior sink

% write Parameters into .txt file if needed.
fileID = fopen([SSParamsPathName, '\Source-Sink-Spiral Detection Parameters.txt'],'w');

fprintf(fileID,'%4s\r\n\r\n','Source-Sink Detection Parameters');
fprintf(fileID,'%20s %15s\r\n\r\n','Parameter','Value');

fprintf(fileID,'%20s %15d\r\n','Nmin',Nmin);
fprintf(fileID,'%20s %15.5f\r\n','Lmin_source',Lmin_source);
fprintf(fileID,'%20s %15.5f\r\n','Lmax_sink',Lmax_sink);
fprintf(fileID,'%20s %15d\r\n','Nnested_source',Nnested_source);
fprintf(fileID,'%20s %15d\r\n','Nnested_sink',Nnested_sink);
fclose(fileID);


