S1=0;


load('D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Wg_Pitot\LC1_01\Pitot_WireWG\matFiles\wireWg.mat')
s=wireWg_m;
b=ones(1,20)/20;
sf=filtfilt(b,1,s);
var_s=filtfilt(b,1,sf.^2);
std_s=sqrt(var_s);
std_s=filtfilt(b,1,sqrt(var_s));
S1=S1+std_s;

load('D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Wg_Pitot\LC1_02\Pitot_WireWG\matFiles\wireWg.mat')
s=wireWg_m;
b=ones(1,20)/20;
sf=filtfilt(b,1,s);
var_s=filtfilt(b,1,sf.^2);
std_s=sqrt(var_s);
std_s=filtfilt(b,1,sqrt(var_s));

S1=S1+std_s;

load('D:\CURRENT_PROJECTS\LC\FabMarcNovDec2014\Data\Wg_Pitot\LC1_03\Pitot_WireWG\matFiles\wireWg.mat')
s=wireWg_m;
b=ones(1,20)/20;
sf=filtfilt(b,1,s);
var_s=filtfilt(b,1,sf.^2);
std_s=sqrt(var_s);
std_s=filtfilt(b,1,sqrt(var_s));

S1=S1+std_s;

S1=S1/3;

S1=decimate(S1,200);



