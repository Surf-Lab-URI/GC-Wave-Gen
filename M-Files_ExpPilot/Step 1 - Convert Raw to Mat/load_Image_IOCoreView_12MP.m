function [IM] = load_Image_IOCoreView_12MP(filename)


imagedim = [4096 3072];
fileoffset = 0;  %file header length in bytes using Dynamic Studio


fid = fopen( filename , 'r' , 'a' );
fseek( fid , fileoffset , 'bof' );
IM = fread( fid , imagedim , 'uint16' , 'a' );
IM=IM';
IM=fliplr(IM); % Flip so flow goes left to right

        
fclose(fid);
%cd(d)