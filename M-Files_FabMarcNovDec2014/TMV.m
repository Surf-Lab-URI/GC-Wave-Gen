clear
clc
% close all

LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
expDirs=dir(LONG); % Find all the experiment directories. 
expDirs=expDirs(3:end-1); %Exclude . and .. dirs as well as the results dir.

ii = 1;

exp_name=expDirs(ii).name;

IRDir = dir([LONG exp_name '/IRMat/']);
IRDir = IRDir(3:end-1);
IRPath = IRDir(1);
load([IRPath.folder,'/',IRPath.name],'IR');
IRSample = IR;
IRfps = 1/(0.1389/6);

usurf = nan(length(IRDir)-1,1);
usurfpx = nan(length(IRDir)-1,1);
usurf0 = nan(length(IRDir)-1,1);
usurf1 = nan(length(IRDir)-1,1);
usurf2 = nan(length(IRDir)-1,1);
TMVTech = nan(length(IRDir)-1,1);

parfor i = 1:length(IRDir)-1
    iA = i %140*6+40;
    iB = iA+1;
    IRPathA = IRDir(iA);
    IRPathB = IRDir(iB);
    
    IRA = load([IRPathA.folder,'/',IRPathA.name],'IR');
    IRA = IRA.IR;
    IRB = load([IRPathB.folder,'/',IRPathB.name],'IR');
    IRB = IRB.IR;
    %%
    % figure(1)
    % hold off
    % imagesc(IRA.img,[19.5,20.25])
    % 
    % figure(2)
    % hold off
    % imagesc(IRB.img,[19.5,20.25])
    
    % tuk = tukeywin(size(IRA.img,1),0.05);
    % Tuk = repmat(tuk,1,size(IRA.img,2));
    
    if mod(iA,24) > 1 && mod(iA,24) < 15
        Tmax = 23;
        Tmin = 19.5;
        IRimAClip = IRA.img;
        IRimAClip(IRimAClip > Tmax) = Tmax;
        IRimBClip = IRB.img;
        IRimBClip(IRimBClip > Tmax) = Tmax;
        IRimAClip(IRimAClip < Tmin) = Tmin;
        IRimBClip(IRimBClip < Tmin) = Tmin;
        TMVTech(i) = 0;
    elseif mod(iA,24) <=1
        Tmax = 0.07;
        Tmin = 0;
        IRimAClip = -ExternalForceImage2D(IRA.img(200:end,:), 0, 1, 0, 3);
        IRimBClip = -ExternalForceImage2D(IRB.img(200:end,:), 0, 1, 0, 3);
        TMVTech(i) = 1;
    else
        Tmax = 0.07;
        Tmin = 0.;
        IRimAClip = -ExternalForceImage2D(IRA.img, 0, 1, 0, 3);
        IRimBClip = -ExternalForceImage2D(IRB.img, 0, 1, 0, 3);
        TMVTech(i) = 2;
    end
    
    
    % IRimAClip(IRimAClip < Tmin) = Tmin;
    % IRimBClip(IRimBClip < Tmin) = Tmin;
    
    % IRimAClip = IRimAClip.*Tuk;
    % IRimBClip = IRimBClip.*Tuk;
    
    % IRimAClip = IRimAClip(200:end,:);
    % IRimBClip = IRimBClip(200:end,:);
    
    h = size(IRimAClip,1);
    % m = abs((1:h)-round(h/2))';
    % 
    % fftA = fft2(IRimAClip);
    % fftB = fft2(IRimBClip);
    % fftCorr = fftB .* conj(fftA);
    % 
    % Xcorr = fftshift(real(ifft2(fftCorr)))./sqrt(sum(sum(IRimAClip.^2)))./sqrt(sum(sum(IRimBClip.^2))); % Cross Correlation
    % Xcorr1D = Xcorr(:,round(size(Xcorr,2)/2));
    % figure(4)
    % hold off
    % plot(Xcorr1D/max(m))
    % % pause
    % % hold on
    % % Xcorr1D = 1./(h-m).*Xcorr1D;
    % % plot(Xcorr1D)
    % [~,Xpky] = max(Xcorr1D) % Find max in the cross correlation
    % dely = Xpky - size(IRimAClip,1)/2-1
    
    corr = zeros(1, round(h/9));
    for m = 1:length(corr)
        A = IRimAClip(1:(end-m),:);
        B = IRimBClip((m+1):end,:);
        A = A - mean(A(:));
        B = B - mean(B(:));
        corr(m) = sum(A.*B, "all") / (norm(A(:)) * norm(B(:)));
    end
    [corrmax, mmax] = max(corr);
    
    % figure(1)
    % hold off
    % imagesc(IRimAClip,[Tmin,Tmax])
    % 
    % figure(2)
    % hold off
    % imagesc(IRimBClip,[Tmin,Tmax])
    % 
    % figure(3)
    % plot(corr)
    % 
    % pause(0.01)
    % 
    % figure(4)
    % hold on
    % plot(i, mmax,'.r')
    
    usurfpx(i) = mmax;
    if TMVTech(i) == 0
        usurf0(i) = mmax;
    elseif TMVTech(i) == 1
        usurf1(i) = mmax;
    else
        usurf2(i) = mmax;
    end
end
usurf = usurfpx*IRfps*IRSample.DX;
usurf0 = usurf0*IRfps*IRSample.DX;
usurf1 = usurf1*IRfps*IRSample.DX;
usurf2 = usurf2*IRfps*IRSample.DX;
%% Plotting Results
figure(5)
hold off
t = (1:length(usurf))/IRfps;
plot(t,usurf,'DisplayName','raw')
hold on
plot(t,usurf0,'.','DisplayName','Tech 0, raw')
plot(t,usurf1,'.','DisplayName','Tech 1, raw')
plot(t,usurf2,'.','DisplayName','Tech 2, raw')
plot(t,movmean(usurf,20),'DisplayName','filtered')
plot(t,movmean(usurf0,20,'omitmissing'),'DisplayName','Tech0, filtered')
% plot(t,movmean(usurf1,20,'omitmissing'),'DisplayName','Tech2, filtered')
% plot(t,movmean(usurf1,20,'omitmissing'),'DisplayName','Tech2, filtered')
legend

%% Saving Results
load_path = [LONG exp_name]
files=dir([load_path '/PIVRaw/PIV/*.mat']);
runResultsDir = [load_path '/Results_Surflab/'];
runResultsfname = [runResultsDir exp_name '_results.mat'];

overwriteRunResults = false;

if ~isfolder(runResultsDir)
    mkdir(runResultsDir)
end

% Create results file for all experiments. If it already exists, save the old
% version as a new file with the creation data in the name. Overwrite the
% original with the new output file.
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

USurf = struct();
USurf.usurf = usurf;
USurf.t = t;

%% Save to run results file
matRun.Usurf = USurf;
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

exps{ii}.USurf = USurf;
matCamp.exps = exps;