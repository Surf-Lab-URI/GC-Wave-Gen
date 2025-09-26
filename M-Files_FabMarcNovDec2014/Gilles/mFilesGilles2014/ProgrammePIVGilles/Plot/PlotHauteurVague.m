clear all
close all

Res = 4.4956e-4;

HV1 = load('E:\ComputedVelocities\20140804\LC1_1\ExpLC1_1_HVague.mat')
% HV1 = load('E:\ComputedVelocities\20140804\LC1_2\ExpLC1_2_HVague.mat')

HV2 = load('E:\ComputedVelocities\20140804\LC2_1\ExpLC2_1_HVague.mat')
% HV2 = load('E:\ComputedVelocities\20140804\LC2_2\ExpLC2_2_HVague.mat')

HV3 = load('E:\ComputedVelocities\20140804\LC3_1\ExpLC3_1_HVague.mat')
HV4 = load('E:\ComputedVelocities\20140804\LC4_3\ExpLC4_2_HVague.mat')


xx=[1:700]/7.2;
inter= [1:570];
smo=1000000
H1=Res*smoothn(HV1.HVague,smo);
H1=H1-H1(1)
H2=Res*smoothn(HV2.HVague,smo);
H2=H2-H2(1)
H3=Res*smoothn(HV3.HVague,smo);
H3=H3-H3(1)
H4=Res*smoothn(HV4.HVague,smo);
H4=H4-H4(1)

xx=[1:700]/7.2-21;
inter = [152:570]; % D'apres les donnees pitots le vent commence a soufler au bout de 21s, soit environ l'image 151/152

figure, plot(xx(inter),H1(inter))
hold on, plot(xx(inter),H2(inter),'r'), plot(xx(inter),H3(inter),'k'), plot(xx(inter),H4(inter),'c')
