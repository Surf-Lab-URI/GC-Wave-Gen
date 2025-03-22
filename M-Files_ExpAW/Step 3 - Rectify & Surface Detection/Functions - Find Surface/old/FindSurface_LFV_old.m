function [imSurf] = FindSurface_LFV_old(img, Step,S,expName)
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
    %%%%
    if isempty(locs)
        surface(i) = nan;
        surface2(i) = nan;
        surface3(i) = nan;
        badFrameBool = 1;
    elseif locs(1) + Step > size(img,1) || locs(1) - Step < 1  % bad frame
        badFrameBool = 1;
    else
        gv2 = abs(Eext2(locs(1)-Step:locs(1)+Step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs(1)-Step-1;
        surface2(i) = locs(1);
        surface3(i) = s_gv2+locs(1)-Step-1;
    end
    
% % %     if i>1101 && i<3800 && ~isempty(locs) %%% modify "i" interval to be in the laser sheet area
% % %         if abs(surface2(i-1)-surface2(i))>Jump_thr 
% % %             [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/4);  %, 'MinPeakProminence', 1e-2), 'npeaks',1);
% % %             [~,Ind] = min(abs(surface2(i-1)-locs));
% % %             surface2(i) = locs(Ind);
% % %             surface3(i) = s_gv2+locs(Ind)-Step-1;
% % %             if abs(surface2(i-1)-surface2(i))>Jump_thr*2
% % %                 L = 1000;
% % %                 cc=polyfit(max(i-L,1):i-1,surface2(max(i-L,1):i-1),4);
% % %                 Segm = polyval(cc,max(i-L,1):i);
% % %                 DiffS = diff(Segm);
% % %                 surface2(i) = surface2(i-1)+DiffS(end);
% % %                 surface3(i) = surface3(i-1)+DiffS(end);
% % %                 MarkInt = [MarkInt i];
% % %                 %surface2(i) = round(Segm(end));
% % %             end
% % %         end
% % %     end
    
    clear gv1 gv2 locs

end

% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;

%% Post-Correction for sudden jumps
% Despiking jumps in Shoaling Waves
surf2_corr = despike_jumps(surface3,Jump_thr,MarkInt,'LFV');
surface3 = surf2_corr;

%% Bad Frame 
Length = 1101:3800;
SS = surface(Length);
SS3 = surface3(Length);
idx=isnan(SS);
fit1 = polyfit(Length(~idx),SS(~idx),50);
fit3 = polyfit(Length(~idx),SS3(~idx),50);
Surf_smooth = polyval(fit1,Length(~idx));
Surf3_smooth = polyval(fit3,Length(~idx));
Diff = abs(medmob2(Surf_smooth'-Surf3_smooth',100));
if max(Diff)>300
    badFrameBool = 1;
else
    badFrameBool = 0;
end

%%
imSurf.surface = surface2; %choose between surface2 and surface3 (surface2 smoother but surface3 better to detect breaking and discontinuities)
imSurf.badFrameBool = badFrameBool;

end
