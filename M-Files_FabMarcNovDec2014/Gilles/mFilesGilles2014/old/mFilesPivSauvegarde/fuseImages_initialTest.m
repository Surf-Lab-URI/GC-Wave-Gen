clear all
close all
load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv1\ExpLC1_dt25ms_1_Piv1_0015_a.mat')
load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv2\ExpLC1_dt25ms_1_Piv2_0015_a.mat')

%% up down matching
nS = 15; %number of pixels to shift
%imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
imgPiv2_shifted = [imgPiv2(nS:end,:);3000*ones(nS-1,2048)];
imgPiv1Piv2 = [imgPiv1 imgPiv2_shifted];
imgPiv1Piv2(536+10:536+35,:) = 4000;
figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([0 3000])


%% left right matching

clear all

load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv1\ExpLC1_dt25ms_1_Piv1_0015_a.mat')
load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv2\ExpLC1_dt25ms_1_Piv2_0015_a.mat')
nS = 15; %number of pixels to shift
%imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
imgPiv2_shifted = [imgPiv2(nS:end,:);3000*ones(nS-1,2048)];
% imgPiv2_shifted = imgPiv2_shifted*2;
imgPiv1Piv2 = [imgPiv1 imgPiv2_shifted];
imgPiv1Piv2(:,1798) = 4000;
imgPiv1Piv2(:,2298) = 4000;
imgPiv1Piv2(:,2048) = 4000;
figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([0 1000]);