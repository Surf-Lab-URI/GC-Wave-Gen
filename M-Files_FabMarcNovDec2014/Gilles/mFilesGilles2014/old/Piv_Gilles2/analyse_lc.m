%%
Struct=compVel
%% Affichage resultats
figure,imagesc(Struct.delx_int.*Struct.MASK)
figure,imagesc(Struct.dely_int.*Struct.MASK)

figure,imagesc(Struct.delxOrb(1:255,1:472).*Struct.MASK), colorbar
figure,imagesc(Struct.delyOrb(1:255,1:472).*Struct.MASK), colorbar

dwhatlevel = Struct.dwhatlevel;
figure, imagesc(dwhatlevel), colorbar

dcor = Struct.dcor;
MASK = Struct.MASK;
dcor1 = 50 * ones(size(dcor)); % matrix of 50s
dcor1(~isnan(MASK)) = dcor(~isnan(MASK)); % matrix of 50s where no velocity calculations, dcor values where calculations were made
figure, imagesc(dcor1), caxis([0 1]), colorbar
numel(find(dcor<0.5))/numel(find(dcor1<40))





%%
Turbx=Struct.delx_int(1:255,1:472)-Struct.delxOrb(1:255,1:472);
Turby=Struct.dely_int(1:255,1:472)-Struct.delyOrb(1:255,1:472);
figure,imagesc(Turbx)
figure,imagesc(Turby)

%% Champs de vitesse
figure, quiver(1:10:472,1:5:255,flipud(Struct.delx_int(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)),flipud(Struct.dely_int(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)))
%% Orbitales
figure, quiver(1:10:472,1:5:255,flipud(Struct.delxOrb(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)),flipud(Struct.delyOrb(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)))
%% Difference
figure, quiver(1:10:472,1:5:255,flipud(Turbx(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)),flipud(Turby(1:5:255,1:10:472).* Struct.MASK(1:5:255,1:10:472)))

