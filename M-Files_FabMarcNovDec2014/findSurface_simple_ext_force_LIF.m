% Function findSurface

%% Author:
% Marc Buckley
%% Last update:
% 08/22/2013
%%
function imSurf = findSurface_simple_ext_force_LIF(IM_a, d1)

PIV=flipud(imgaussian(medfilt2(IM_a,[5 5]),5));
clear surfacePIV
for i=1:size(PIV,2)
gv1 = abs(PIV(:,i));
[~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
if (isempty(locs))
    locs=min(find(gv1==max(gv1)));
end
surfacePIV(i) =locs;
end
b=ones(1,32)/32;
surfacePIV=2048-filtfilt(b,1,surfacePIV);
badFrameBool=0;



img=imgaussian(d1,5);
clear surface
for i=1:size(img,2)
gv1 = abs(img(:,i));
[~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
if (isempty(locs))
    locs=find(gv1==max(gv1));
end
surface(i) =locs;
end
b=ones(1,32)/32;
surface=filtfilt(b,1,surface);
badFrameBool=0;
surfTechnique=1;
surface=surface(714:714+2047)-1835+364;
%
   
clear surface2
Wedge=2;
Wline=0;
Wterm=0;
Sigma1=32;
step=1;
Sigma2=step;
%
Eext1 = ExternalForceImage2D(img,Wline, Wedge, Wterm,Sigma1);
Eext2 = ExternalForceImage2D(img,Wline, Wedge, Wterm,Sigma2);
surface2 = nan(size(img,2),1);
%
badFrameBool = 0;
for i=1:size(img,2)
%     i = 10034;
 %i
    gv1 = abs(Eext1(:,i));
    [~,locs] = findpeaks(gv1, 'minpeakheight', max(gv1)/2, 'npeaks',1);
    if isempty(locs)
        surface2(i) = nan;
    elseif locs + step > size(img,1) || locs - step < 1  % bad frame
        badFrameBool = 1;
        break
    else
        gv2 = abs(Eext2(locs-step:locs+step,i));
        [~, s_gv2] = max(gv2);
        surface2(i) = s_gv2+locs-step-1;
    end
    clear gv1 gv2 locs
end

b=ones(1,32)/32;
surface2=filtfilt(b,1,surface2);
surface2=surface2(714:714+2047)'-1835+364;
surfTechnique2=2;

if (abs(mean(surfacePIV-surface))<abs(mean(surface-surface2)))
imSurf.surface = surface;
imSurf.badFrameBool = badFrameBool;
imSurf.surfTechnique = surfTechnique;
else
imSurf.surface = surface2;
imSurf.badFrameBool = badFrameBool;
imSurf.surfTechnique = surfTechnique2;
end

end