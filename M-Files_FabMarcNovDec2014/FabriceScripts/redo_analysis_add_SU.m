% TREAT VELOCITY
clear

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
%Load PIV Surf
fprintf(['pair ' num2str(image_pair_number)])

load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
%Surface detection and Creating Masks
d1=imresize(imgPivsurf,176.77/103.48); %Resizing to match PIV
imSurf = findSurface_simple_ext_force((medfilt2(d1)), 1);
Surface_PIV=imSurf.surface;
Surface_PIV=Surface_PIV(1:end-1); %just so that it's even number for FFT
fprintf(['.'])
%Load PIV data
load([load_path '\PIVMat\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.mat']); %replace ~ with path

%Calculated Orbital Vel
pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;

%pivRes.pf_surf = Surface_PIV(X) -1838+370;
%altitude=[0 compVel.zPIV-compVel.GS/2]; %compVel.GS/2 on average with moving surface
%X=724+compVel.GS-1:compVel.GS:724+2047-compVel.GS;
fprintf(['.'])
transfo = generateTransfo_LC_noLFV( compVel, Surface_PIV, pivRes); % 0 is there to compare SU(0,:) with surface
%
SU = transfo.SU; 
SU = SU -1838+370;
pivRes.pf_surf =SU(1,:); %THE SURFACE EXACTELY
SU = SU(2:end,:); % all but surface; % first line is zeta=0 THE SURFACE EXACTELY

ORBX = transfo.ORBX;
ORBX = ORBX(2:end,:);
ORBZ = transfo.ORBZ;
ORBZ = ORBZ(2:end,:);
%
pivRes.GS = compVel.GS;
pivRes.zPIV = compVel.zPIV;



u = compVel.delx.*compVel.mask;
w = compVel.dely.*compVel.mask;
%
intrp_u = transformVelField_decay_forFab( u, pivRes, SU );
intrp_w = transformVelField_decay_forFab( w, pivRes, SU );
%% Subtract ORBs
intU_minusORBX = intrp_u - ORBX;
intW_minusORBZ = intrp_w - ORBZ;
%%
%% inverse transform
% pivRes.zPIV = 2:4:2047;
u_turb = reverseTransformVelField_decay_forFab( intU_minusORBX, pivRes, SU );
w_turb = reverseTransformVelField_decay_forFab( intW_minusORBZ, pivRes, SU );
uuTest = reverseTransformVelField_decay_forFab( intrp_u, pivRes, SU);

intrp_u_turb = transformVelField_decay_forFab( u_turb, pivRes, SU );
intrp_w_turb = transformVelField_decay_forFab( w_turb, pivRes, SU );

compVel.u=u;
compVel.w=w;
compVel.u_turb=u_turb;
compVel.w_turb=w_turb;
compVel.intrp_u=intrp_u;
compVel.intrp_w=intrp_w;
compVel.ORBX=ORBX;
compVel.ORBZ=ORBZ;
compVel.SU=SU;
compVel.pf_surf=pivRes.pf_surf;
compVel.intrp_u_turb=intrp_u_turb;
compVel.intrp_w_turb=intrp_w_turb;
%field = ['delx'; 'dely']; compVel = rmfield(compVel,field);
field = ['INTdelx'; 'INTdelz']; compVel = rmfield(compVel,field);


outfile = [load_path '\PIVMat_TURB\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)];
save(outfile, 'compVel', 'pivRes' );disp([' done.']);
fprintf('\r')



end


end


%FROM compVel Getting a transformed and reverse transform 
pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;
pivRes.pf_surf=compVel.pf_surf;
%d=reverseTransformVelField_decay_forFab(compVel.ORBX, pivRes, compVel.SU )


