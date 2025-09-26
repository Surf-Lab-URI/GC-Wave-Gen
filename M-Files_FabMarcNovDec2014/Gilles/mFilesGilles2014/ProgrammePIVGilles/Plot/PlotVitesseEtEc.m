clear all
close all


A1 = load('\\Beo\data\ExpLC2_dt7ms_1\ProcessedVelocities\ExpLC2_dt7ms_1_calcEcOrb2.mat');
A1Bis = load('\\Beo\data\ExpLC2_dt7ms_3\ProcessedVelocities\ExpLC2_dt7ms_3_calcEcOrb.mat');
A3 = load('\\Beo\data\ExpLC3_dt7ms_1\ProcessedVelocities\ExpLC3_dt7ms_1_calcEcOrb2.mat');
A4 = load('\\Beo\data\ExpLC4_dt7ms_1\ProcessedVelocities\ExpLC4_dt7ms_1_calcEcOrb2.mat');
PivRes = 40d-6;
DeltaT= 7e-3;
Sc = PivRes/DeltaT;

% xx = ([1:300]+200)/7.2;
xx = ([1:300]+200)/7.2-21;

yy = [1:150]*8*PivRes;

figure, imagesc(xx,yy,(A1.EcAndEcOrb.delx_mean(1:150,1:300)*Sc).^2+(A1.EcAndEcOrb.dely_mean(1:150,1:300)*Sc).^2), colorbar



figure, plot(A1.EcAndEcOrb.delx_mean(1:150,70)*Sc,-yy)
hold on, plot(A1.EcAndEcOrb.delx_mean(1:150,100)*Sc,-yy, 'c'), plot(A1.EcAndEcOrb.delx_mean(1:150,150)*Sc,-yy,'m'), plot(A1.EcAndEcOrb.delx_mean(1:150,182)*Sc,-yy,'--r'),  plot(A1.EcAndEcOrb.delx_mean(1:150,211)*Sc,-yy,'--k')













figure, imagesc(xx,yy,(A1.EcAndEcOrb.delxTurb_mean(1:150,1:300)*Sc).^2+(A1.EcAndEcOrb.delyTurb_mean(1:150,1:300)*Sc).^2), colorbar





figure, imagesc(A1.EcAndEcOrb.ec), colorbar
figure, imagesc((A1.EcAndEcOrb.delxTurb_mean.^2+A1.EcAndEcOrb.delyTurb_mean.^2)/2)
figure, imagesc(A3.EcAndEcOrb.ec)
figure, imagesc(A4.EcAndEcOrb.ec)

figure, surf(smoothn(A1.EcAndEcOrb.ec))
figure, surf((A1.EcAndEcOrb.delxTurb_mean.^2+A1.EcAndEcOrb.delyTurb_mean.^2)/2)

figure, surf(smoothn(A3.EcAndEcOrb.ec))
figure, surf(smoothn(A4.EcAndEcOrb.delx_mean))






