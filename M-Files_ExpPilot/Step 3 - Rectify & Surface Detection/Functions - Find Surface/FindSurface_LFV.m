function [imSurf] = FindSurface_LFV(img, Step,S,expName)
% This function finds the surface in PIV Surf as well as LFV images. It was
% basically written by Marc Buckley on the August 22 2013, and the original
% file can be found in the following directory:  \MFiles\surface_detection\
% functions.
%
%    This function implements the basic snake segmentation contour. A snake
%    is an active (moving) contour, in which the points are being attracted
%    by edges and other image boundaries. In order to keep  contour smooth,
%    a membrane and thin plate energy are used as contour regularization.

%
% outputs,
%  Eextern : The energy function described by the image
%

Sigma1 = 20;   % Sigma used to calculated image derivatives
Sigma2 = Step; % Sigma used to calculated image derivatives

Eext1 = ExternalForceImage2D_fab(img,Sigma1).*S;
Eext2 = ExternalForceImage2D_fab(img,Sigma2);
Eext1(Eext1<0)=0;
badFrameBool = 0;
warning off
Jump_thr = 15; % error if surface jumps by more than <Jump_thr> pixels
HighWindExp = [7:12,19:24,31:36,43:48];
if ismember(str2double(expName),HighWindExp)
    Jump_thr = 15*2;
end
MarkInt = []; % marked interval to be re-interpolated
for i=1:size(img,2)
    
    gv1 = abs(Eext1(:,i));
    
    %     [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4, 'npeaks',1);
    
    %%%% Modified by Fabio A. for Shoaling waves
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2); % ,max(gv1)/2.5, 'npeaks',1);
    if i >2000 && i<2500 %%% if rotated PIVSurf
        [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4); % ,max(gv1)/2, 'npeaks',1);
    end
    
    if str2double(expName) > 36
        [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4); % ,max(gv1)/2, 'npeaks',1);
    end
    
    if str2double(expName) > 48
        [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/6); % ,max(gv1)/2, 'npeaks',1);
    end
    %%%%
    if isempty(locs)
        surface(i) = nan;
        surface2(i) = nan;
        badFrameBool = 1;
    elseif locs(1) + Step > size(img,1) || locs(1) - Step < 1  % bad frame
        badFrameBool = 1;
    else
        gv2 = abs(Eext2(locs(1)-Step:locs(1)+Step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs(1)-Step-1;
        surface2(i) = locs(1);
    end
    
    clear gv1 gv2 locs
    
end

% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;

%%  Remove Points with no light
if str2double(expName)<49
    idx2 = isnan(surface');
    x_s = (1:length(surface));
    s_s = surface;
    s_s(idx2) = 1;
    LinInd = sub2ind([size(img,1),size(img,2)],s_s',x_s');
    LinInd2 = find(img(round(LinInd))<1);
    LinInd2 = unique([LinInd2,LinInd2+50,LinInd2-50]);
    LinInd2(LinInd2<1) = [];
    LinInd2(LinInd2>length(surface)) = [];
    surface(LinInd2) = NaN;
    idx3 = isnan(surface);
    L(1) = find(idx3==0,1);
    L(2) = length(idx3)-find(fliplr(idx3)==0,1);
    surface(1:L(1)-1) = surface(L(1));
    surface(L(end)+1:end) = surface(L(end));
    idx4 = isnan(surface);
    surface = interp1(x_s(~idx4),surface(~idx4),x_s,'pchip','extrap');
end
%% Post-Correction for sudden jumps
% Despiking jumps in Shoaling Waves
surf2_corr = despike_jumps(surface,Jump_thr,'LFV',expName);
surface3 = surf2_corr;

%% Last attempt to recover frame before marking as bad frame
if str2double(expName)<49
    
    Length = 801:3900;
    
    if str2double(expName) > 36
        Surf = surface3;
        Surf(901:1200) = NaN;
        Surf = Surf(Length);
        Idx = isnan(Surf);
        fitS = polyfit(Length(~Idx),Surf(~Idx),10);
        surface3(901:1200) = polyval(fitS,Length(101:400));
    end
    
    Range = [min(setdiff(1:length(surface),LinInd2)) max(setdiff(1:length(surface),LinInd2))];
    surface3([1:Range(1)-1,Range(end)+1:end]) = NaN;
    
    % Interpolate data with no lights
    if Range(end)<Length(end)
        Ind_fit2 = max(setdiff(1:length(surface),LinInd2))-500:max(setdiff(1:length(surface),LinInd2));
        fit2 = polyfit(Ind_fit2,surface3(Ind_fit2),2);
        Tail = smooth(polyval(fit2,Ind_fit2(end):Length(end)),length(Ind_fit2)/2);
        surface3(Ind_fit2(end)+1:Length(end)) = Tail(2:end)+surface3(Ind_fit2(end))-Tail(1);
    end
    if Range(1)>Length(1)
        Ind_fit2 = min(setdiff(1:length(surface),LinInd2)):min(setdiff(1:length(surface),LinInd2))+500;
        fit2 = polyfit(Ind_fit2,surface3(Ind_fit2),2);
        Tail = smooth(polyval(fit2,Length(1):Ind_fit2(1)),length(Ind_fit2)/2);
        surface3(Length(1):Ind_fit2(1)-1) = Tail(1:end-1)+surface3(Ind_fit2(1))-Tail(end);
        
    end
    
    SS = surface(Length);
    [~,locs1] = findpeaks(abs(diff(SS)), 'minpeakheight', Jump_thr);locs1(locs1<Length(1) | locs1>Length(end)) = [];
    SS = surface3(Length);
    for ii = 1:length(locs1)
        SS(locs1(ii)-70:locs1(ii)+70) = NaN;
    end
    idx = isnan(SS);
    fit1 = polyfit(Length(~idx),SS(~idx),20);
    Surf_smooth = polyval(fit1,Length);
    Surf_smooth = [surface3(1:Length(1)+200-1) Surf_smooth(201:end-200) surface3(Length(end)-200+1:end)];
    Diff = abs(Surf_smooth-surface3);
    Int = 1:length(Diff);
    Int(Diff<Jump_thr) = [];
    Int(isnan(surface3(Int))) = [];
    IntP = Int+100;
    IntM = Int-100;
    Int = unique([Int,IntP,IntM]);
    %
    for i = 1:length(Int)
        gv1 = abs(Eext1(:,Int(i)));
        [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/6, 'MinPeakProminence', 1e-2); %, 'npeaks',1);
        if ~isempty(locs)
            [~,I] = min(abs(locs-Surf_smooth(Int(i))));
            gv2 = abs(Eext2(locs(I)-Step:locs(I)+Step,Int(i)));
            [~, s_gv2] = max(gv2);
            surface3(Int(i)) = s_gv2+locs(I)-Step-1;
        end
    end
    surf2_corr = despike_jumps(surface3,Jump_thr,'LFV',expName);
    surface3 = surf2_corr;
end

%% Bad Frame
Length = 801:3900;
if str2double(expName)>48 && str2double(expName)<55
    Length = 401:3600;
end
SS = surface(Length);
[~,locs1] = findpeaks(abs(diff(SS)), 'minpeakheight', Jump_thr);
SS = surface3(Length);
for ii = 1:length(locs1)
    SS(max(1,locs1(ii)-100):min(length(SS),locs1(ii)+100)) = NaN;
end
idx = isnan(SS);
L2(1) = find(idx==0,1);
L2(2) = length(idx)-find(fliplr(idx)==0,1);
SS(1:L2(1)-1) = surface(Length(1)+(1:L2(1)-1));
SS(L2(end)+1:end) = surface(Length(1)+(L2(end)):end-((length(surface)-Length(end))));
idx2 = isnan(SS);
Length2 = 1101:3500;
fit1 = polyfit(Length(~idx2),SS(~idx2),20);
Surf_smooth = polyval(fit1,Length2);

Diff = abs(Surf_smooth'-surface3(Length2)');
[~,Imax] = max(Diff);
if max(Diff)>4*Jump_thr && Imax < 1000 && Imax > 3400
    badFrameBool = 1;
else
    badFrameBool = 0;
end

%%
imSurf.surface = surface3; %choose between surface2 and surface3 (surface2 smoother but surface3 better to detect breaking and discontinuities)
imSurf.badFrameBool = badFrameBool;

end
