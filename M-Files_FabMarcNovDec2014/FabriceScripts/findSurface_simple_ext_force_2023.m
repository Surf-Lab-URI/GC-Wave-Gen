% Function findSurface

%% Author:
% Marc Buckley
%% Last update:
% 08/22/2013
%%
function imSurf = findSurface_simple_ext_force(img, step)
%
Options.Wedge=2;
Options.Wline=0;
Options.Wterm=0;
Options.Sigma1=50;
Options.Sigma2=step;
%
Eext1 = ExternalForceImage2D(img,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma1);
Eext2 = ExternalForceImage2D(img,Options.Wline, Options.Wedge, Options.Wterm,Options.Sigma2);
surface = nan(size(img,2),1);
%
badFrameBool = 0;
for i=1:size(img,2)
%     i = 10034;

 warning off
    gv1 = abs(Eext1(:,i));
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
    if isempty(locs)
        surface(i) = nan;
    elseif locs + step > size(img,1) || locs - step < 1  % bad frame
        badFrameBool = 1;
        break
    else
        gv2 = abs(Eext2(locs-step:locs+step,i));
        [~, s_gv2] = max(gv2);
        surface(i) = s_gv2+locs-step-1;
    end
    clear gv1 gv2 locs
end
inn = ~isnan(surface);
i1 = (1:numel(surface)).';
pp = interp1(i1(inn),surface(inn),'linear','pp');
surface_s = fnval(pp,linspace(i1(1),i1(end),length(surface)))';

b=ones(1,32)/32;
surface=filtfilt(b,1,surface_s);
%
% imSurf = surface;
% surface(abs(surface-nanmean(surface))>4*std(surface)) = nan;
imSurf.surface = surface;
imSurf.badFrameBool = badFrameBool;
end