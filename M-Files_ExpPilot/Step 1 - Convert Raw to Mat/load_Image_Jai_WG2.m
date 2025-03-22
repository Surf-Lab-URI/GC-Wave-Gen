function [img] = load_Image_Jai_WG2(imagename)

nl = 59;
fileoffset = 0;  %file header length in bytes using Dynamic Studio
fid = fopen( imagename);
fseek( fid , fileoffset , 'bof' );
img = fread(fid,nl*1600,'uint16');
img = reshape(img,1600,nl);
        
fclose(fid);
%cd(d)