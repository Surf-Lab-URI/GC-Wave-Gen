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

% Results arrays
usurfpx    = nan(length(IRDir)-1,1);  % Selected technique's result (px/frame)
usurfpx0   = nan(length(IRDir)-1,1);  % TMVTech 0 result for every frame (px/frame)
usurfpx1   = nan(length(IRDir)-1,1);  % TMVTech 1 result for every frame (px/frame)
TMVTech    = nan(length(IRDir)-1,1);  % Which technique the gate selected (0 or 1)
corrmax0   = nan(length(IRDir)-1,1);  % Peak correlation for TMVTech 0 (weighted avg)
corrmax1   = nan(length(IRDir)-1,1);  % Peak correlation for TMVTech 1
Ndots_used = nan(length(IRDir)-1,1);  % Number of dots that passed threshold

% Parameters
Ndots = 8;
dot_spacing_px = 47.29;
x_dot1 = 44;
T_threshold = 20.75;
corr_threshold = 0.5;  % Minimum peak correlation to accept a dot's measurement

% Change parfor to for when debugging with figures
parfor i = 1:900%410%1:1200
    iA = i;
    iB = iA+1;
    IRPathA = IRDir(iA);
    IRPathB = IRDir(iB);
    
    IRA = load([IRPathA.folder,'/',IRPathA.name],'IR');
    IRA = IRA.IR;
    IRB = load([IRPathB.folder,'/',IRPathB.name],'IR');
    IRB = IRB.IR;
    
    % figure(4)
    % hold off
    % imagesc(IRA.img,[19.5,20.25])
    % 
    % figure(5)
    % hold off
    % imagesc(IRB.img,[19.5,20.25])
    
    % =====================================================================
    % TMVTech 0: dot-based cross-correlation on clipped temperature images
    % =====================================================================
    Tmax = 23;
    Tmin = 19.5;
    IRimAClip = IRA.img;
    IRimAClip(IRimAClip > Tmax) = Tmax;
    IRimBClip = IRB.img;
    IRimBClip(IRimBClip > Tmax) = Tmax;
    IRimAClip(IRimAClip < Tmin) = Tmin;
    IRimBClip(IRimBClip < Tmin) = Tmin;
    
    h = size(IRimAClip,1);
    
    weighted_sum_mmax = 0;
    sum_weights = 0;
    n_dots_ok = 0;
    
    for n = 1:Ndots
        col_start = round(x_dot1 + (n-1)*dot_spacing_px);
        col_end   = round(x_dot1 + n*dot_spacing_px);
        IRimAStrip = IRimAClip(:, col_start:col_end);
        IRimBStrip = IRimBClip(:, col_start:col_end);
        
        if max(IRimAStrip,[],'all') > T_threshold && max(IRimBStrip,[],'all') > T_threshold
            Nm = round(h/12) + 1;             % include m=0
            corr = zeros(1, Nm);
            m_vals = 0:(Nm-1);
            for k = 1:Nm
                m = m_vals(k);
                if m == 0
                    A = IRimAStrip;
                    B = IRimBStrip;
                else
                    A = IRimAStrip(1:(end-m),:);
                    B = IRimBStrip((m+1):end,:);
                end
                A = A - mean(A(:));
                B = B - mean(B(:));
                corr(k) = sum(A.*B, "all") / (norm(A(:)) * norm(B(:)));
            end
            [cmax, kmax] = max(corr);
            
            if cmax > corr_threshold
                kmax_sub = subpixel_peak(corr, kmax);
                mmax_sub = kmax_sub - 1;      % convert index -> shift
                weighted_sum_mmax = weighted_sum_mmax + cmax * mmax_sub;
                sum_weights = sum_weights + cmax;
                n_dots_ok = n_dots_ok + 1;
                
                % figure(1)
                % hold off
                % imagesc(IRimAStrip,[Tmin,Tmax])
                % daspect([1,1,1])
                % 
                % figure(2)
                % hold off
                % imagesc(IRimBStrip,[Tmin,Tmax])
                % daspect([1,1,1])
                % 
                % figure(3)
                % plot(m_vals, corr)
                % 
                % mmax_sub
            end
        end
    end
    
    if sum_weights > 0
        dots_present = true;

        mmax0 = weighted_sum_mmax / sum_weights;
        cmax0 = sum_weights / n_dots_ok;

        % Plot displaced image pair for debugging
        % figure(4)
        % hold off
        % imagesc(IRA.img(1:(end-round(mmax0)),:),[19.5,21])
        % 
        % figure(5)
        % hold off
        % imagesc(IRB.img((1+round(mmax0)):end,:),[19.5,21])
        % pause
    else
        dots_present = false;
        mmax0 = nan;
        cmax0 = nan;
    end


    % =====================================================================
    % TMVTech 1: filtered full-image cross-correlation (always compute)
    % =====================================================================
    IRimAFilt = -ExternalForceImage2D(IRA.img(1:end,:), 0, 1, 0, 6);
    IRimBFilt = -ExternalForceImage2D(IRB.img(1:end,:), 0, 1, 0, 6);
    
    h1 = size(IRimAFilt,1);
    Nm = round(h1/12) + 1;
    corr = zeros(1, Nm);
    m_vals = 0:(Nm-1);
    for k = 1:Nm
        m = m_vals(k);
        if m == 0
            A = IRimAFilt;
            B = IRimBFilt;
        else
            A = IRimAFilt(1:(end-m),:);
            B = IRimBFilt((m+1):end,:);
        end
        A = A - mean(A(:));
        B = B - mean(B(:));
        corr(k) = sum(A.*B, "all") / (norm(A(:)) * norm(B(:)));
    end
    [cmax1, kmax] = max(corr);
    kmax_sub = subpixel_peak(corr, kmax);
    mmax1 = kmax_sub - 1;                     % convert index -> shift
    
    % figure(1)
    % hold off
    % imagesc(IRimAFilt)
    % daspect([1,1,1])
    % 
    % figure(2)
    % hold off
    % imagesc(IRimBFilt)
    % daspect([1,1,1])
    % 
    % figure(3)
    % plot(m_vals, corr)
    
    % =====================================================================
    % Store both results plus the gate decision
    % =====================================================================
    usurfpx0(i)   = mmax0;
    usurfpx1(i)   = mmax1;
    corrmax0(i)   = cmax0;
    corrmax1(i)   = cmax1;
    Ndots_used(i) = n_dots_ok;
    
    if dots_present
        TMVTech(i) = 0;
        usurfpx(i) = mmax0;
    else
        TMVTech(i) = 1;
        usurfpx(i) = mmax1;
    end
end

% Convert from px/frame to physical velocity units
usurf  = usurfpx  * IRfps * IRSample.DX;
usurf0 = usurfpx0 * IRfps * IRSample.DX;
usurf1 = usurfpx1 * IRfps * IRSample.DX;


function m_sub = subpixel_peak(corr, k_int)
% Parabolic fit through three points around the peak for sub-pixel accuracy.
% Returns the sub-pixel-refined *index* into corr (still 1-based); the
% caller is responsible for converting index -> shift if needed.
    if k_int <= 1 || k_int >= length(corr)
        m_sub = k_int;
        return
    end
    y1 = corr(k_int - 1);
    y2 = corr(k_int);
    y3 = corr(k_int + 1);
    denom = (y1 - 2*y2 + y3);
    if denom == 0
        m_sub = k_int;
    else
        delta = 0.5 * (y1 - y3) / denom;
        m_sub = k_int + delta;
    end
end

%%
% Example comparison plot (run after the parfor loop finishes)
figure
hold on
plot(usurf0, '.-', 'DisplayName', 'TMVTech 0 (dots)')
plot(usurf1, '.-', 'DisplayName', 'TMVTech 1 (filtered)')
xlabel('Frame index')
ylabel('Surface velocity (m/s)')
legend

figure
plot(usurf0, usurf1, '.')
hold on
plot(xlim, xlim, 'k--')  % 1:1 line
xlabel('TMVTech 0 velocity (m/s)')
ylabel('TMVTech 1 velocity (m/s)')
axis equal
%% Plotting Results
figure(5)
hold off
t = (1:length(usurf))/IRfps;
plot(t,usurf,'DisplayName','raw')
hold on
plot(t,usurf0,'.','DisplayName','Tech 0, raw')
plot(t,usurf1,'.','DisplayName','Tech 1, raw')
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