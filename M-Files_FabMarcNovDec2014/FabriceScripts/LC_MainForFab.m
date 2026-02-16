 %Main for LC using Fab input
%Marc 09/24/2015
%% load Fab input
load('C:\Users\Asilab\Desktop\Downloads\fab_LC_Im.mat');
%% generate decaying surface and orbital matrixes on full PS length (to avoid edge effects on PIV part)
% my function doesn't like uneven surface vectors, hence the "end-1"
% altitude = 0:4:20;
%100 sec on my machine
pivRes.zPIV = compVel.zPIV;
pivRes.xPIV = compVel.xPIV;
pivRes.GS = compVel.GS;
pivRes.mask = compVel.mask;

%pivRes.pf_surf = Surface_PIV(X) -1838+370;
%altitude=[0 compVel.zPIV-compVel.GS/2]; %compVel.GS/2 on average with moving surface
%X=724+compVel.GS-1:compVel.GS:724+2047-compVel.GS;

tic, transfo = generateTransfo_LC_noLFV( compVel, Surface_PIV, pivRes); toc, % 0 is there to compare SU(0,:) with surface
%

SU = transfo.SU; 
SU = SU(2:end,:); % all but surface;
SU = SU -1838+370;
ORBX = transfo.ORBX;
ORBX = ORBX(2:end,:);
ORBZ = transfo.ORBZ;
ORBZ = ORBZ(2:end,:);
%
pivRes.GS = compVel.GS;
pivRes.zPIV = compVel.zPIV;
pivRes.pf_surf =SU(1,:);

usmth = smoothn(compVel.delta_x, 0, 'robust');
u = usmth.*compVel.mask;
wsmth = smoothn(compVel.delta_z, 0, 'robust');
w = wsmth.*compVel.mask;
%
intrp.u = transformVelField_decay_forFab( u, pivRes, SU );
intrp.w = transformVelField_decay_forFab( w, pivRes, SU );
%% Subtract ORBs
intU_minusORBX = intrp.u - ORBX;
intW_minusORBZ = intrp.w - ORBZ;
%%
%% inverse transform
% pivRes.zPIV = 2:4:2047;
uu = reverseTransformVelField_decay_forFab( intU_minusORBX, pivRes, SU );
ww = reverseTransformVelField_decay_forFab( intW_minusORBZ, pivRes, SU );
uuTest = reverseTransformVelField_decay_forFab( intrp.u, pivRes, SU);
% 
