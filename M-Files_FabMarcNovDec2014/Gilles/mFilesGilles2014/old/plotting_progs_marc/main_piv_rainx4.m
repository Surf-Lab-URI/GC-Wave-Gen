clear all
close all

imfolder = '\\Afsx1\piv2\Exp3';
dt = 200e-6;
pixres = 47.4e-6;
cd(imfolder);
piv_velocity('Exp3_', pixres, dt)
