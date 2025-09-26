clear all
close all
pixres = 47.4e-6;

%% user:
imfolder = '\\Afsx1\PIV6\Exp23\';
savefolder = '\\Afsx1\PIV_processed\'
files_prefix = 'Exp23_';
dt = 130e-6;
startnum = 1;
endnum = 1000;
% direc = imfolder;
%%

% cd(imfolder);
piv_velocity(files_prefix, pixres, dt,startnum,endnum, imfolder, savefolder)
