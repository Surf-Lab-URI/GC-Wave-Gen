clear
clc

LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';


expDirs=dir(LONG); % Find all the experiment directories. 
expDirs=expDirs(3:end-1); %Exclude . and .. dirs as well as the results dir.
i = 1;
% for i = 1:length(expDirs)
ii = 2; %Index of the experiment run
num_of_digits = 3; %Number of digits used in image number in file names
exp_name = expDirs(ii).name;
load_path = [LONG exp_name];
files=dir([load_path '/PIVRaw/PIV/*.mat']);
runResultsDir = [load_path '/Results_Surflab/'];
runResultsfname = [runResultsDir exp_name '_results.mat'];
number_of_pair=length(files)/2;

overwriteRunResults = false;

if ~isfolder(runResultsDir)
    mkdir(runResultsDir)
end

resultsfname = [LONG, 'Results_Surflab/results.mat'];
if ~isfile(resultsfname)
    save(resultsfname,'-v7.3')
else
    dirResultsFName = dir(resultsfname);
    ResultsfnameDate = datetime(dirResultsFName.date);
    ResultsfnameDate.Format = 'yyyy-MM-dd_HHmmss';

    copyfile(resultsfname,[LONG, 'Results_Surflab/results_' char(ResultsfnameDate) '.mat'])
end
matCamp = matfile(resultsfname,'Writable',true);

% Check if run results file exists. If it doesn't, create it. If it does,
% open it but save a copy first just in case something bad happens.
if ~isfile(runResultsfname) || overwriteRunResults
    save(runResultsfname,'-v7.3');
else
    dirRunResultsFName = dir(runResultsfname);
    runResultsfnameDate = datetime(dirRunResultsFName.date);
    runResultsfnameDate.Format = 'yyyy-MM-dd_HHmmss';

    copyfile(runResultsfname,[runResultsDir exp_name '_results_' char(runResultsfnameDate) '.mat'])
end

matRun = matfile(runResultsfname,'Writable', true);
%%
pairs = cell(number_of_pair,1);

tic
for image_pair_number = 0:number_of_pair-1
    pair = struct();
    aPath = [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat'];
    bPath = [load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat'];
    pair.imSurfa = FindSurfaceCapillary(aPath);
    pair.imSurfa.path = aPath;
    pair.imSurfb = FindSurfaceCapillary(bPath);
    pair.imSurfb.path = bPath;
    pairs{image_pair_number+1} = pair;
    disp(['Found surface for pair ' num2str(image_pair_number)])
end
toc
%% Build struct with surfaces and timing
tic

dx = 1/17697.69; %m per pix
pps = 7.2; % pairs per second
spp = 1/pps; % seconds per pair
dt_pair = 10e-3; %delay between images in a pair
t = nan(number_of_pair*2,1);
surfs = nan(number_of_pair*2,length(pairs{1}.imSurfa.surface));
for image_pair_number = 0:number_of_pair-1
    t(image_pair_number*2 + 1) = 2*image_pair_number*spp;
    t(image_pair_number*2 + 2) = 2*image_pair_number*spp+dt_pair;

    if ~isempty(pairs{image_pair_number+1})
        surfs(image_pair_number*2 + 1,:) = pairs{image_pair_number+1}.imSurfa.surface;
        surfs(image_pair_number*2 + 2,:) = pairs{image_pair_number+1}.imSurfb.surface;
        fprintf('wrote image pair %d to surfs matrix\n',image_pair_number)
    else
        fprintf('No surface for image_pair_number = %d\n',image_pair_number)
    end
end
eta = (surfs-mean(surfs(1:20,:),'all'))*dx;
x_eta = (0:(size(eta,1)-1))*dx;

toc
%% Make Struct to save surfs
Surfs = struct();
Surfs.surfs = surfs;
Surfs.dx = dx;
Surfs.t = t;
Surfs.pps = pps;
Surfs.spp = spp;
Surfs.dt_pair = dt_pair;
Surfs.eta = eta;
Surfs.x_eta = x_eta;
%% Save surfs in run results file
matRun.Surfs = Surfs;
%% Load previously saved surfs
Surfs = matRun.Surfs;
%% Save surfs to campaign results file
matCampVars = who(matCamp);
exps = cell(length(expDirs),1);
exps{ii} = struct();
if ~ismember('exps',matCampVars)
    matCamp.exps = cell(length(expDirs),1);
elseif ~iscell(matCamp.exps)
    matCamp.exps = cell(length(expDirs),1);
else
    exps = matCamp.exps;
end

exps{ii}.Surfs = Surfs;
matCamp.exps = exps;
% end