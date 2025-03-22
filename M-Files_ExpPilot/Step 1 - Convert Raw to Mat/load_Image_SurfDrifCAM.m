function [img] = load_Image_SurfDrifCAM(imagename)

nl = 2048;
fileoffset = 0;  %file header length in bytes using Dynamic Studio
fid = fopen( imagename);
fseek( fid , fileoffset , 'bof' );
img = fread(fid,nl*2048,'uint16');
img = reshape(img,2048,nl);
        
fclose(fid);
%cd(d)