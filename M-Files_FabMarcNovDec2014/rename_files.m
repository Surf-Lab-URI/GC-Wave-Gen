% TREAT VELOCITY

LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\']
DIRS=dir(LONG);
DIRS=DIRS(3:end);

for ii=1:length(DIRS)

exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = ['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\' exp_name];
write_path =['H:\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\' exp_name];
files=dir([load_path '\PIVMat_TURB\*.mat']);
number_of_im=length(files);

for image_pair_number=1:number_of_im;
%Load PIV Surf
fprintf(['pair ' num2str(image_pair_number)])

load([load_path '\PIVMat_TURB\' files(image_pair_number).name]); %replace ~ with path

outfile = [write_path '\PIVMat_TURB\' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number-1)];
save(outfile, 'compVel', 'pivRes' );
disp([' done.']);
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


