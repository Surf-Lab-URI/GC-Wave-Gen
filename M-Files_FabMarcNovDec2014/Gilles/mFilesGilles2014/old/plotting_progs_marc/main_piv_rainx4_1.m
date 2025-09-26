clear all
close all
pixres = 47.4e-6;

%% user:
imfolder = '\\Afsx1\PIV6\Movie13';
exp = 'Exp25_';
dt = 130e-6;
startnum = 1;
endnum = 1000;
%%

cd(imfolder);
piv_velocity(exp, pixres, dt,startnum,endnum)
