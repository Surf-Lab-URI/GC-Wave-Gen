% Fonction fuseImages
%% Object:
% Corrige et fusionne deux images de PIV
% Correction des differences d'angle, de gain et d'offset entre les cameras
% Fusion des deux images par gradient lineaire
% Ici, le vent souffle de droite a gauche
%% Arguments: 
% Image1, Image2, tform 
%% Resultat: 
% experiment number, image pair number, image letter, fused image
%% Auteurs:
% Gilles Bouille, Fabrice Veron, Marc Buckley
%% Derniere mise a jour:
% 14/08/2013
%% Exemple:

% u = [ 65 246 2048]'; v = [ 617 1786 582]'; x = [1795 1991 3785]'; y = [632 1808 594]'; tform = maketform('affine',[u v],[x y]);
% fusedIm = fuseImagesLC(imgPiv1,imgPiv2,tform);
% ans = 
%        fused_im: [2040x3785 double]
%        exp_name: '3'
%     im_pair_num: 440
%       im_letter: 'a'

%%
function PivFuse_cut = fuseImages_lc(imgPiv1,imgPiv2,tform)
%% Correction due a l'angle et placement relatif des deux images
imgPiv1_t=imtransform(imgPiv1,tform,'Xdata',[1 3785],'Ydata',[1 2048], 'FillValues',-1); 
imgPiv2_big=zeros(size(imgPiv1_t))-1; imgPiv2_big(1:2048,1:2048)=imgPiv2;

%% Correction de l'offset et gain 
imgPiv2_corr = imgPiv2*1+144;

%% Fusion des deux images par gradient lineaire sur la zone de recouvrement
PIV1=imgPiv1_t*0;
PIV1(imgPiv1_t==-1)=1;
PIV2=imgPiv2_big*0;
PIV2(imgPiv2_big==-1)=1;
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
    PivFuse(i,j)=imgPiv2_corr(i,j);
    end
    for j=position(i)+1:2048
         PivFuse(i,j)=((j-position(i))/recouvrement(i))*imgPiv1_t(i,j)+((2048-j)/recouvrement(i))*imgPiv2_corr(i,j);
    end
    for j=2049:3785
        PivFuse(i,j)=imgPiv1_t(i,j);
    end
end
%% Rognure des bords pour obtenir une image rectngulaire
PivFuse_cut=PivFuse(:,:);
