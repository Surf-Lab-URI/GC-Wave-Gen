function [IM1] = load_Image_IOCoreView_48MP(filename)


imagedim = [7920 6004];
fileoffset = 0;  %file header length in bytes using Dynamic Studio


fid = fopen( filename , 'r' , 'a' );
fseek( fid , fileoffset , 'bof' );
IM1 = fread( fid , imagedim , 'uint16' , 'a' );
IM1=IM1';
        
fclose(fid);
%cd(d)