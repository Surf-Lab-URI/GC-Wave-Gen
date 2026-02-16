LONG=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz']
DIRS=dir(LONG);
DIRS=DIRS(3:end);

PITOT=['D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Wg_Pitot']
DIRS_PITOT=dir([PITOT '\LC*'])


ii=1;  %GO FROM 1 TO 12 ONLY!!!!
exp_name=DIRS(ii).name;
exp_name_Pitot=DIRS_PITOT(ii).name;

num_of_digits = 3;
load_path = [LONG,'\' exp_name];
files=dir([load_path '\PIVRaw\PIV\*.mat']);
number_of_pair=length(files)/2;

%constants
DX=5.65571e-05; %Resolution for PIV image in m/pix

L=[];    

for image_pair_number=100:200
 image_pair_number
%PIV LFV to get the surface
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
imgPivsurf1 = imgPivsurf;
load([load_path '\PIVRaw\PIVSURF\' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
imgPivsurf2 = imgPivsurf;

%Detecting the surface to calculate wavenumber
d1=imresize(imgPivsurf1,176.77/103.48); %Resizing to match PIV
imSurf1 = findSurface_simple_ext_force((medfilt2(d1)), 1); 
surface_profile=imSurf1.surface*DX; surface_x=[1:length(imSurf1.surface)]*DX; %scaling  surface profile in m
%d2=imresize(imgPivsurf2,176.77/103.48);
%imSurf2 = findSurface_simple_ext_force((medfilt2(d2)), 1);
wave_numb=wav_num(surface_x,surface_profile,DX);   %calculation of wavenumber from surface profile
wave_amplitude=std(surface_profile);


s=12;
 
load([PITOT,'\' exp_name_Pitot, '\Pitot_WireWG\matFiles\pitot.mat' ]); %replace ~ with path
Time_of_image_pair=round((50+s+(image_pair_number)/432*60)*200);
U10=pit_m_s(Time_of_image_pair-100:Time_of_image_pair+100);
U10=mean(U10);
        
La=Langmuir_number(wave_amplitude,wave_numb,U10);

L=[L La];
end

        
        




