LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\']
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1:length(DIRS)

exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [LONG exp_name];
files=dir([load_path '\PIVRaw\PIV\*.mat']);
number_of_pair=length(files)/2;

for image_pair_number=0:number_of_pair-1
%PIV
load([load_path '\PIVRaw\PIV\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imgPiv;
load([load_path '\PIVRaw\PIV\' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = imgPiv;
 
%PIV Surf
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurf1 = imgPivsurf;
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurf2 = imgPivsurf;

%Surface detection and Creating Masks
U1 = [147 49;2024 57; 1995 1004; 161 999];
X1 = [147 49; 2024 49; 2024 1004; 147 1004];
T1 = fitgeotrans(U1,X1,'projective');
d1=imwarp(imgPivsurf1,T1,'cubic');
d2=imwarp(imgPivsurf2,T1,'cubic');

d1=imresize(d1,176.9769/105.5880);%Resizing to match PIV
d2=imresize(d2,176.9769/105.5880);%Resizing to match PIV
s1=d1(30:3525,755:755+2047); %croping
s2=d2(30:3525,755:755+2047);

imSurf1 = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
imSurf2 = findSurface_simple_ext_force_2023((medfilt2(s2)), 1);
imSurf1.surface=imSurf1.surface-1716+287;
imSurf2.surface=imSurf2.surface-1716+287;
[h, w] = size(IM_a); %image height and width
mask1=ones(size(IM_a));
mask2=ones(size(IM_a));
for i=1:w
    mask1(1:round(imSurf1.surface(i)),i)=NaN;
    mask2(1:round(imSurf2.surface(i)),i)=NaN;
end
warning off
imSurf1.mask=mask1;
imSurf2.mask=mask2;

IntrWndw=[256 128 64 32 16 8];
GrdSpc=[128 64 32 16 8 4];
compVel =  computeVelocities_marc_quick_nofilt(IM_a, IM_b, mask1, mask2, IntrWndw, GrdSpc);
compVel.DX=1/17697.69; %m per pix
compVel.DT=10d-3; % sec per image pair -  DELTA_T= 10 milisec

outfile = [load_path '\PIVMat_2023\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'compVel', 'imSurf1', 'imSurf2');
disp(['pair ' num2str(image_pair_number) ' velocity done.']);

%%
imSurf = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
Surface_PIV=imSurf.surface;

pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;

transfo = generateTransfo_LC_noLFV_2023( compVel, Surface_PIV, pivRes); % 0 is there to compare SU(0,:) with surface

SU = transfo.SU; 
SU = SU(2:end,:); % all but surface; % first line is zeta=0 THE SURFACE EXACTELY
SU = SU -1716+287;
ORBX = transfo.ORBX;
ORBX = ORBX(2:end,:);
ORBZ = transfo.ORBZ;
ORBZ = ORBZ(2:end,:);
%
pivRes.GS = compVel.GS;
pivRes.zPIV = compVel.zPIV;
pivRes.pf_surf =SU(1,:);

u = compVel.delx.*compVel.mask;
w = compVel.dely.*compVel.mask;
%
intrp_u = transformVelField_decay_forFab( u, pivRes, SU );
intrp_w = transformVelField_decay_forFab( w, pivRes, SU );

intU_minusORBX = intrp_u - ORBX;
intW_minusORBZ = intrp_w - ORBZ;

u_turb = reverseTransformVelField_decay_forFab( intU_minusORBX, pivRes, SU );
w_turb = reverseTransformVelField_decay_forFab( intW_minusORBZ, pivRes, SU );
uuTest = reverseTransformVelField_decay_forFab( intrp_u, pivRes, SU);
