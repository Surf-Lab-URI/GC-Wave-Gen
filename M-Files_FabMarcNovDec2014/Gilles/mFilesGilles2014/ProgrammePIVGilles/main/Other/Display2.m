clear all
close all
num_of_digits=4;
image_pair_number =364;

%% Chargement compVel
load(['E:\ComputedVelocities\20140730\LC1dt10\ExpLC1dt10_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
LC1dt10 = compVel;
load(['E:\ComputedVelocities\20140730\LC1dt15\ExpLC1dt15_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
LC1dt15 = compVel;
load(['E:\ComputedVelocities\20140730\LC3dt10\ExpLC3dt10_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
LC3dt10 = compVel;
load(['E:\ComputedVelocities\20140730\LC3dt15\ExpLC3dt15_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat'])
LC3dt15 = compVel;
clear compVel;

%% Deplacements transversaux
figure, imagesc(LC1dt10.delx_ints.*LC1dt10.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC1dt15.delx_ints.*LC1dt15.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC3dt10.delx_ints.*LC3dt10.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC3dt15.delx_ints.*LC3dt15.MASK), colorbar, caxis([-5 5])

%% Deplacements verticaux

figure, imagesc(LC1dt10.dely_ints.*LC1dt10.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC1dt15.dely_ints.*LC1dt15.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC3dt10.dely_ints.*LC3dt10.MASK), colorbar, caxis([-5 5])
figure, imagesc(LC3dt15.dely_ints.*LC3dt15.MASK), colorbar, caxis([-5 5])

%% Correlation

figure, imagesc(LC1dt10.dcor.*LC1dt10.MASK), colorbar
figure, imagesc(LC1dt15.dcor.*LC1dt15.MASK), colorbar
figure, imagesc(LC3dt10.dcor.*LC3dt10.MASK), colorbar
figure, imagesc(LC3dt15.dcor.*LC3dt15.MASK), colorbar


dcor = LC1dt10.dcor; MASK = LC1dt10.MASK; dcor1 = 50 * ones(size(dcor)); dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); numel(find(dcor<0.5))/numel(find(dcor1<40))
dcor = LC1dt15.dcor; MASK = LC1dt15.MASK; dcor1 = 50 * ones(size(dcor)); dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); numel(find(dcor<0.5))/numel(find(dcor1<40))
dcor = LC3dt10.dcor; MASK = LC3dt10.MASK; dcor1 = 50 * ones(size(dcor)); dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); numel(find(dcor<0.5))/numel(find(dcor1<40))
dcor = LC3dt15.dcor; MASK = LC3dt15.MASK; dcor1 = 50 * ones(size(dcor)); dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); numel(find(dcor<0.5))/numel(find(dcor1<40))

%% Niveau de calcul

figure, imagesc(LC1dt10.dwhatlevel.*LC1dt10.MASK), colorbar
figure, imagesc(LC1dt15.dwhatlevel.*LC1dt15.MASK), colorbar
figure, imagesc(LC3dt10.dwhatlevel.*LC3dt10.MASK), colorbar
figure, imagesc(LC3dt15.dwhatlevel.*LC3dt15.MASK), colorbar

nanmean(nanmean(LC1dt10.dwhatlevel.*LC1dt10.MASK))
nanmean(nanmean(LC1dt15.dwhatlevel.*LC1dt15.MASK))
nanmean(nanmean(LC3dt10.dwhatlevel.*LC3dt10.MASK))
nanmean(nanmean(LC3dt15.dwhatlevel.*LC3dt15.MASK))