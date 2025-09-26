clear all
close all

%% Ouverture des fichiers
% Image calib
matFramePiv1_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\piv1\Movie11_Scene8_piv1_0.raw',2048,2048);
matFramePiv2_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\piv2\Movie11_Scene8_piv2_0.raw',2048,2048);
matFramePivsurf_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\pivsurf\Movie11_Scene8_pivsurf_0.raw',2048,2048);
matFrameLfv_Sc8 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene8\lfv\Movie11_Scene8_lfv_0.raw',2048,2048);

imgPiv1=matFramePiv1_Sc8.img;
imgPiv2=matFramePiv2_Sc8.img;
imgPivsurf=flipud(fliplr(matFramePivsurf_Sc8.img));

% u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [564 559 1562]'; tform = maketform('affine',[u v],[x y]);
% imgPivsurf_t=imtransform(imgPivsurf,tform,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-10);

% Exemples PIV
matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_0001.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_0001.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_0001.raw',2048,2048);

matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_2001.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_2001.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_2001.raw',2048,2048);

matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_2002.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_2002.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_2002.raw',2048,2048);

matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_2003.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_2003.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_2003.raw',2048,2048);

matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_2004.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_2004.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_2004.raw',2048,2048);

matFramePiv1_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv1\Movie11_Scene20_piv1_2005.raw',2048,2048);
matFramePiv2_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\piv2\Movie11_Scene20_piv2_2005.raw',2048,2048);
matFramePivsurf_Sc20 = saveToMatSingleFrame('F:\data\rawPIV\Movie11_Scene20\pivsurf\Movie11_Scene20_pivsurf_2005.raw',2048,2048);


imgPiv1=matFramePiv1_Sc20.img;
imgPiv2=matFramePiv2_Sc20.img;
imgPivsurf=flipud(fliplr(matFramePivsurf_Sc20.img));

%% Fusion et correction Piv1 Piv2
u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
imgPiv1_t=imtransform(imgPiv1,tform,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1); %mettre -0.1 puis faire un arrondi a l'unite pour ne pas modifier les  valeurs
B=zeros(size(imgPiv1_t))-1; B(1:2048,1:2048)=imgPiv2*2; C=B+imgPiv1_t;
% figure, imagesc(C), colormap(bone), caxis([0 500])
% %  Controle 
% figure, imagesc(B), colormap(bone), caxis([0 500]), axis([1600 2200 800 1200])
% figure, imagesc(imgPiv1_t), colormap(bone), caxis([0 500]), axis([1600 2200 800 1200])
% 
% figure, imagesc(B), colormap(bone), caxis([0 500]), axis([1600 2200 1800 2048])
% figure, imagesc(imgPiv1_t), colormap(bone), caxis([0 500]), axis([1600 2200 1800 2048])

% % Correction Offset et gain 
imgPiv2_144 = imgPiv2+144;
% [N1,X1] = hist(imgPiv1(:),1000);
% [N2,X2] = hist(imgPiv2_144(:),1000);
% figure, plot(X1,N1,'k', X2, N2, 'r')

% Fusion de l'image avec gradient 
PIV1=imgPiv1_t*0;
PIV1(imgPiv1_t==-1)=1;
PIV2=B*0;
PIV2(B==-1)=1;
PIVF= PIV1+PIV2;

for i=1:2048
   position(i)=0;
   for j=1:2047
       position(i)=position(i)+PIVF(i,j);
   end
   recouvrement(i)=2048-position(i)-1;
end
PivFuse=imgPiv1_t*0;
for i=1:2048
    for j=1:position(i)
    PivFuse(i,j)=imgPiv2_144(i,j);
    end
    for j=position(i)+1:2048
         PivFuse(i,j)=((j-position(i))/recouvrement(i))*imgPiv1_t(i,j)+((2048-j)/recouvrement(i))*imgPiv2_144(i,j);
    end
    for j=2049:3785
        PivFuse(i,j)=imgPiv1_t(i,j);
    end
end
% figure, imagesc(PivFuse), colormap(bone), caxis([0 500]), colorbar


%% Correction et superposition Pivsurf, PivFuse 

%Pas parfais mais adjuge 
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
%Donnees initiales
% u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [544 559 1562]'; tformSurf = maketform('affine',[u v],[x y]);
% imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);
% figure, imagesc(imgPivsurf_t), colormap(bone)

%% Correction et superposition Lfv, PivFuse 

%% Espace travail

%Pas parfais mais adjuge 
u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1756]'; y = [552 554 1562]'; tformSurf = maketform('affine',[u v],[x y]);
imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);


% u = [ 598 1577 990]'; v = [ 1089 1090 1479]'; x = [744 3262 1762]'; y = [549 553 1559]'; tformSurf = maketform('affine',[u v],[x y]);
% imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);

% % figure, imagesc(imgPivsurf_t)
% 
% 
% u =[403 412 1577 1575]'; v=[1068 1455 1076 1490]'; x=[244 241 3262 3274]'; y=[566 1599 560 1577]';
% tformSurf =cp2tform([u v],[x y],'projective');
% imgPivsurf_t=imtransform(imgPivsurf,tformSurf,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1);


% figure, imagesc(imgPivsurf), colormap(bone)
% 
% 
% figure, imagesc(PivFuse), colormap(bone), caxis([0 500])
% figure, imagesc(imgPivsurf_t), colormap(bone)
% 
% figure, imagesc(PivFuse), colormap(bone), caxis([0 800]), axis([500 1500 400 800])
% figure, imagesc(imgPivsurf_t), axis([500 1500 400 800])
% 
% 
% figure, imagesc(PivFuse), colormap(bone), caxis([0 800]), axis([3500 3785 400 800])
% figure, imagesc(imgPivsurf_t), axis([3500 3785 400 800])





img = imgPivsurf_t;
smth_vert = nan(size(img));
grad_vert = nan(size(img));
mask = nan(size(img));
surface = nan(1,size(img,2));
%
for i=1:size(img,2)
%     keyboard
    imgi = img(:,i);
    %locate outliers and nan them
    imgistd = std(imgi);
    imgimean = mean(imgi);
    imgi(abs(imgi-imgimean)>3*imgistd)=nan;
    %smooth each column
    smth_vert(:,i) = smoothn(imgi,10000);
    %compute gradient on each column
    grad_vert(:,i) = gradient(smth_vert(:,i));
    gv = grad_vert(:,i);
    [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/2, 'npeaks',1);
    surface(i) = locs;
    clear imgi gv locs pks
    clear grad_vert %%%%%%%%%%%%%%%%%%%%%%%%
end
%
%locate outliers and nan them
surfstd = std(surface);
surfmean = mean(surface);
surface(abs(surface-surfmean)>3*surfstd)=nan;
smth_surface = smoothn(surface);


smth_surface = smoothn(surface,1000);




figure, imagesc(PivFuse), colormap(bone), caxis([0 800])
% figure, imagesc(imgPivsurf_t), colormap(bone)
hold on , plot(smth_surface)



















decoupage(1:600,:)=1*imgPivsurf_t(1:600,:);
decoupage(601:2048,:)=10*PivFuse(601:2048,:);
figure, imagesc(decoupage), colormap(bone)
















































%% Old


% clear all
% close all
% load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv1\ExpLC1_dt25ms_1_Piv1_0015_a.mat')
% load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv2\ExpLC1_dt25ms_1_Piv2_0015_a.mat')
% 
% %% up down matching
% nSv = 15; %number of pixels to shift vertical
% 
% %imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
% imgPiv2_shifted = [imgPiv2(nSv:end,:);3000*ones(nSv-1,2048)];
% imgPiv1Piv2 = [imgPiv2_shifted imgPiv1 ];
% 
% figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([0 1000])
% 
% 
% %% left right matching
% 
% clear all
% 
% load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv1\ExpLC1_dt25ms_1_Piv1_0015_a.mat')
% load('F:\data\ExpLC1_dt25ms_1\rawImages\Piv2\ExpLC1_dt25ms_1_Piv2_0015_a.mat')
% nS = 15; %number of pixels to shift
% %imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
% imgPiv2_shifted = [imgPiv2(nS:end,:);3000*ones(nS-1,2048)];
% % imgPiv2_shifted = imgPiv2_shifted*2;
% imgPiv1Piv2 = [imgPiv1 imgPiv2_shifted];
% imgPiv1Piv2(:,1798) = 4000;
% imgPiv1Piv2(:,2298) = 4000;
% imgPiv1Piv2(:,2048) = 4000;
% figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([0 1000]);
% 
% %%
% 
% 
% 
% nSv = 15; %number of pixels to shift vertical
% nSh =290;
% %imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
% imgPiv2_shifted = [imgPiv2(nSv:end,:);3000*ones(nSv-1,2048)];
% imgPiv1_shifted = [imgPiv1(:,nSh:2048)];
% imgPiv1Piv2 = [imgPiv2_shifted imgPiv1_shifted ];
% 
% figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([0 300])
% 

% %imgPiv2_shifted = [nan(nS-1,2048); imgPiv2(nS:end,:)];
% imgPiv2_shifted = [imgPiv2(nSv:end,:);3000*ones(nSv-1,2048)];
% imgPiv1_shifted = [imgPiv1_resized_vert(:,nSh:2048)];
% imgPiv1Piv2_1 = [2*imgPiv2_shifted imgPiv1_shifted ];
% 
% % figure, imagesc(imgPiv1Piv2), colormap(bone), caxis([100 500])
% figure, imagesc(imgPiv1Piv2_1), colormap(bone), caxis([0 500])
% 
% figure, imagesc(imgPiv1_t),colormap(bone), caxis([0 500])
%   
% %la mem image en crouant de l'autre cote
% 
% imgPiv2_shifted = [imgPiv2(nSv:end,:);3000*ones(nSv-1,2048)];
% imgPiv2_shifted2 = imgPiv2_shifted(:,1:2048-nSh+1);
% imgPiv1Piv2_2 = [2*imgPiv2_shifted2 imgPiv1_resized_vert ];
% figure, imagesc(imgPiv1Piv2_2), colormap(bone), caxis([0 500])
% 

% %test
% 
% IM1=imgPiv1Piv2_1(800:2048-nSv+1,1733:2049);
% IM2=imgPiv1Piv2_2(800:2048-nSv+1,1733:2049);
% mask=ones(2034,317);
% compVel =  computeVelocities('toto', 0, IM1, IM2, mask);
% 
% figure, quiver(1:7:78,1:20:309,compVel.delta_x(1:20:309,1:7:78),compVel.delta_z(1:20:309,1:7:78))
% 
% 
% figure, imagesc(compVel.delta_x), colorbar
% figure, imagesc(compVel.delta_z), colorbar
% 

