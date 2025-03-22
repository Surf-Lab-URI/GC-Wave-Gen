function [BadFramePIVSurf,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface,PIVSurf_PIVMatch, index] = ExtractSurface_PIVSurf_Rob(norm_ptv_surf,PIVSurf_CRR,fused_image_rotate_crop)


Y1 = abs(round(PIVSurf_CRR.YPos(1))); Y2 = Y1 + size(fused_image_rotate_crop,1) - 1;
X1 = abs(round(PIVSurf_CRR.XPos(1))); X2 = X1 + size(fused_image_rotate_crop,2) - 1;
%Coordinate matching PIV image

h = fspecial('gaussian',16,10);
%PIVSurf_CRR.img(:,4645:end)=PIVSurf_CRR.img(:,4645:end)-11; %offset correction for second tap (right side of image)
PIVSurf_PIVMatch=PIVSurf_CRR.img(:,X1:X2);
% NS=norm_ptv_surf(:,X1:X2);
% f=NS(1,:); NS=[repmat(f,14095,1)];
% PIVSurf_PIVMatch=PIVSurf_PIVMatch./NS;
PIVSurf_PIVMatch=imfilter(PIVSurf_PIVMatch,h,'replicate');
%PIVSurf image exactely matching PIV images
imshow(PIVSurf_PIVMatch, [])

%PIVSurf_Strip = PIVSurf_CRR.img(Y1+1500:Y2, 1+361:end-361);  %Correct size for matlab 2014
%%%ATTENTION%%% use the following instead if using Matlab version 2015 and
%%%above! NORM_PIV_SURF was obtained with a version 2014 and is larger by
%%%two columns - There was a change in "IMROTATE" routine use line 34 of
%%%CorrectPIVSurfv_v2.
%%NN=361 for old Matlab NN=360 for new one
 
%imshow(PIVSurf_PIVMatch, [])

% BW = PIVSurf_PIVMatch;
% 
% BW(PIVSurf_PIVMatch <80) = 0;
% BW(BW>1) = 1;
% imshow(BW, [])
% 
% EDGE =gradient(BW);
% for o = 1:length(EDGE(:,1))
%     Surf(o) = nansum(abs(EDGE(o,:))); 
% end
% index = nanmin(find(Surf == nanmax(Surf)))

% NN=360;
% PIVSurf_Strip = PIVSurf_CRR.img(Y1+1500:Y2, 1+NN:end-NN);  %Correct size for matlabe 2015 and later
% PIVSurf_Strip=PIVSurf_Strip./NORM_PIV_SURF;
% PIVSurf_Strip=imfilter(PIVSurf_Strip,h,'replicate');

% PIVSurf image exactely matching height of PIV images, but longer
% 361 gets rid of NaN from the rotation and matching of PIVSurf with PIV images
%1500 is a vertical offset to get closer to the surface and reduce image
%size in which to detect the surface

%PF_XPos_in_PIVSurf = X1-360:X2-360; %to re-place surface from strip into PIV coordinates

%Creating filter around surface
%PIVSurf_PIVMatch=imgaussian(PIVSurf_PIVMatch,200);
index = 7500;
PIVSurf_Strip = PIVSurf_PIVMatch(nanmin(nanmax(index-4050, 51),length(PIVSurf_PIVMatch(:,1)) -8050) :nanmin(nanmax(index-4050, 51)+8000,length(PIVSurf_PIVMatch(:,1))-50) , :);
PIVSurf_Strip_T = PIVSurf_PIVMatch(nanmin(nanmax(index-4000, 101),length(PIVSurf_PIVMatch(:,1)) -8000) :nanmin(nanmax(index-4000, 101)+8000,length(PIVSurf_PIVMatch(:,1))) , :);
PIVSurf_Strip_B = PIVSurf_PIVMatch(nanmin(nanmax(index-4100, 1),length(PIVSurf_PIVMatch(:,1)) -8100) :nanmin(nanmax(index-4100, 1)+8000,length(PIVSurf_PIVMatch(:,1))-100) , :);


%PIVSurf_Strip = PIVSurf_PIVMatch(nanmin(nanmax(index-3050, 51),length(PIVSurf_PIVMatch(:,1)) -6050) :nanmin(nanmax(index-3050, 51)+6000,length(PIVSurf_PIVMatch(:,1))-50) , :);
% PIVSurf_Strip_T2 = PIVSurf_PIVMatch(nanmin(nanmax(index-2900, 301),length(PIVSurf_PIVMatch(:,1)) -6000) :nanmin(nanmax(index-2900, 301)+6000,length(PIVSurf_PIVMatch(:,1))) , :);
% PIVSurf_Strip_B2 = PIVSurf_PIVMatch(nanmin(nanmax(index-3200, 1),length(PIVSurf_PIVMatch(:,1)) -6100) :nanmin(nanmax(index-3200, 1)+6000,length(PIVSurf_PIVMatch(:,1))-100) , :);
% PIVSurf_Strip_T=PIVSurf_Strip_T./NS(nanmin(nanmax(index-3050, 51),length(PIVSurf_PIVMatch(:,1)) -6050) :nanmin(nanmax(index-3050, 51)+6000,length(PIVSurf_PIVMatch(:,1))-50), :);
% PIVSurf_Strip_B=PIVSurf_Strip_B./NS(nanmin(nanmax(index-3050, 51),length(PIVSurf_PIVMatch(:,1)) -6050) :nanmin(nanmax(index-3050, 51)+6000,length(PIVSurf_PIVMatch(:,1))-50), :);
S_T=(imfilter(PIVSurf_Strip_T,fspecial('gaussian',64,16),'replicate'));
S_B=(imfilter(PIVSurf_Strip_B,fspecial('gaussian',64,16),'replicate'));





% S_T2=(imfilter(PIVSurf_Strip_T2,fspecial('gaussian',64,16),'replicate'));
% S_B2 = (imfilter(PIVSurf_Strip_B2,fspecial('gaussian',64,16),'replicate'));
% 
% 
% 
% %%
% PIVSurf_Strip_T3 = PIVSurf_PIVMatch(nanmin(nanmax(index-2800, 501),length(PIVSurf_PIVMatch(:,1)) -6000) :nanmin(nanmax(index-2800, 501)+6000,length(PIVSurf_PIVMatch(:,1))) , :);
% PIVSurf_Strip_B3 = PIVSurf_PIVMatch(nanmin(nanmax(index-3300, 1),length(PIVSurf_PIVMatch(:,1)) -6100) :nanmin(nanmax(index-3300, 1)+6000,length(PIVSurf_PIVMatch(:,1))-100) , :);
% % PIVSurf_Strip_T=PIVSurf_Strip_T./NS(nanmin(nanmax(index-3050, 51),length(PIVSurf_PIVMatch(:,1)) -6050) :nanmin(nanmax(index-3050, 51)+6000,length(PIVSurf_PIVMatch(:,1))-50), :);
% % PIVSurf_Strip_B=PIVSurf_Strip_B./NS(nanmin(nanmax(index-3050, 51),length(PIVSurf_PIVMatch(:,1)) -6050) :nanmin(nanmax(index-3050, 51)+6000,length(PIVSurf_PIVMatch(:,1))-50), :);
% S_T3=(imfilter(PIVSurf_Strip_T3,fspecial('gaussian',64,16),'replicate'));
% S_B3=(imfilter(PIVSurf_Strip_B3,fspecial('gaussian',64,16),'replicate'));
%%
S=((S_T-S_B));
% S2= ((S_T2-S_B2));
% S3= ((S_T3-S_B3));
% S4 = S2.*S.*S3;
imshow(S, [])

S_positive = S;
S_negative = S;
S_positive(S_positive<0) = 0;
S_negative(S_negative>0) = 0;

% imshow(S_negative, [])
% imshow(S_positive, [])
S_check = [zeros(100,3931); S_negative(100:end-1, :) + S_positive(1:end-100, :)];

% S_check = S_check+[zeros(100,3931); S_negative(99:end-2, :) + S_positive(1:end-100, :)];
% S_check = S_check+[zeros(101,3931); S_negative(10:end-1, :) + S_positive(1:end-101, :)];
% S_check = S_check+[zeros(100,3931); S_negative(100:end-1, :) + S_positive(1:end-100, :)];
figure(3)
imshow(S_check, [])
SE = strel('disk',2);
S3 =  imdilate(S_check,SE);
imshow(S3, [])
S34 = imbinarize(S3, 9);
SE2 = strel('sphere',20);
S35 =  imdilate(S34,SE2);

%figure
[BW5, L, N, A] = bwboundaries(S35);  %10 GOOD
index10 = find(cellfun('length', BW5)>4000)

bw_6 = zeros(size(S_check));
for i = 1:length(index10)
for h = 1:nanmax(cellfun('length', BW5))
    X_check(h) = BW5{index10(i),1}(h,2); 
    Y_check(h) = BW5{index10(i),1}(h,1);
    bw_6(BW5{index10(i),1}(h,2),BW5{index10(i),1}(h,1)) = 1;
end
end
f = fit(X_check', Y_check', 'poly2')
% param.minMajorAxis = 1;
% param.maxMajorAxis = 10000;
% param.rotation = 0;
% param.rotationSpan = 180;
% param.minAspectRatio =  .0001;
% param.numBest = 1; 
% param.randomize = 3; 
%     
% bestFits_guess_surf = ellipseDetection(bw_6,param);
coeff_fit = coeffvalues(f);
 for i = 1:3931
     guess_line(i) =coeff_fit(3) + coeff_fit(2)*i + coeff_fit(1).*(i^(2));
 end
% for h = 1:length(BW5)
%    if (nanmin(abs(h - index2)) == 0)
%   plot(BW5{h,1}(:,2),BW5{h,1}(:,1), '.')
%   hold on
%    else
%    end
% end
%  
% plot(guess_line, '.k')
clear index2
for g = 1:length(guess_line)
    for h = 1:length(BW5)
       dis_part(h) = nanmin(sqrt((g-BW5{h,1}(:,2)).^(2) + (guess_line(g) - BW5{h,1}(:,1)).^(2)));
    end
    index2(g) = find(dis_part ==nanmin(dis_part));
end

% index2 = find(cellfun('length', BW5) == nanmax(cellfun('length', BW5)));

%index2 = [index2; find(cellfun('length', BW5) > 300)];
D = size(S)
S_BW = zeros(D(1), D(2));
for L = 1:length(BW5)
    if (nanmin(abs(L - index2)) == 0)
        for k = 1:cellfun('length', BW5(L))
            S_BW(BW5{L,1}(k,1), BW5{L,1}(k,2)) = 1;
        end
    end 
end
S_BW2 = imfill(S_BW);
 S2 = S_check.*S_BW2;
%  imshow(S2, [])
%  SE = strel('disk',40);
% S3 =  imdilate(S2,SE);
% imshow(S3, [])
PIVSurf_Surface_Tem = FindSurface2_rob(S2,5,S2); %5
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.')
% pause(5)
[fitresult, gof] = nyqist_guess_rob(find(isnan(PIVSurf_Surface_Tem.surface) == 0), PIVSurf_Surface_Tem.surface(isnan(PIVSurf_Surface_Tem.surface) == 0));
if (gof.rsquare >0.35)
     params =  coeffvalues(fitresult);
     phase_longwave_ptv_surf = params(3);
     wavenum_longwave = params(2);
     offset = params(4);
     wave_height = params(1);
     SURF_TEST = PIVSurf_Surface_Tem.surface;
     SURF_TEST(isnan(PIVSurf_Surface_Tem.surface) == 0) = wave_height*cos(wavenum_longwave*PIVSurf_Surface_Tem.surface(isnan(PIVSurf_Surface_Tem.surface) == 0) +  phase_longwave_ptv_surf) + offset;
     PIVSurf_Surface_Tem.surface = SURF_TEST;
else
%finds surface in the strip
% imshow(S, [])
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.')
index_NaN = (find(isnan(PIVSurf_Surface_Tem.surface))); 
gradient_check = nanmean(diff(PIVSurf_Surface_Tem.surface));
gradient2_check = nanmedian(diff(diff(PIVSurf_Surface_Tem.surface))); 
%gradient3_check = nanmean(diff(diff(diff(PIVSurf_Surface_Tem.surface)))); 
index_NaN = (find(isnan(PIVSurf_Surface_Tem.surface)));
SURF_TEST = PIVSurf_Surface_Tem.surface;
try
for j = 1:length(index_NaN)
   
    try
        SURF_TEST(index_NaN(j)) = gradient_check*(index_NaN(j)-(nanmax(index_NaN+1)))+ (gradient2_check./2).*((index_NaN(j)- (nanmax(index_NaN+1)))^2) + SURF_TEST((nanmax(index_NaN+1)));  %+ (gradient2_check./2).*((index_NaN(j)- (nanmax(index_NaN+1)))^2) +  (gradient3_check./6).*((index_NaN(j)- (nanmax(index_NaN+1)))^3)
    catch
        SURF_TEST(index_NaN(j)) = gradient_check*(index_NaN(j)-(nanmin(index_NaN-1))) + (gradient2_check./2).*((index_NaN(j)- (nanmin(index_NaN-1)))^2) + SURF_TEST((nanmin(index_NaN-1))); %(gradient2_check./2).*((index_NaN(j)- (nanmin(index_NaN-1)))^2) + (gradient3_check./6).*((index_NaN(j)- (nanmin(index_NaN-1)))^3)+ 
    end  
end 
PIVSurf_Surface_Tem.surface = SURF_TEST;
catch
end
end
% try
% % imshow(S, [])
% hold on 
% plot(PIVSurf_Surface_Tem.surface-50, '.')
%  catch
%  end
BadFramePIVSurf=PIVSurf_Surface_Tem.badFrameBool;

PIVSurf_Surface_Raw = PIVSurf_Surface_Tem.surface+nanmin(nanmax(index-4050, 51),length(PIVSurf_PIVMatch(:,1) -8050)) -50 ;

PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
%plot(PIVSurf_Surface_Raw-(nanmin(nanmax(index-4050, 51),length(PIVSurf_PIVMatch(:,1) -8050)) -50), '.')
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
%plot(PIVSurf_Surface_Raw-(nanmin(nanmax(index-4050, 51),length(PIVSurf_PIVMatch(:,1) -8050)) -50), '.')
% try
%    % PIVSurf_Surface_Raw = filt_drop_rob(PIVSurf_Surface_Raw, PIVSurf_Strip,nanmin(nanmax(index-4050, 51),length(PIVSurf_PIVMatch(:,1) -8050)));
% 
% catch
%     warning('on');
%     warning(['Problem using FiltSpray function. Ignoring the ' ...
%              'FiltSpray function for PIVSurf image']);
%     PIVSurf_Surface_Raw = PIVSurf_Surface_Raw;
%     BadFramePIVSurf=1;
%  end
%         
% try
%     PIVSurf_Surface_Raw = filt_spray(PIVSurf_Surface_Raw);
% 
% catch
%     warning('on');
%     warning(['Problem using FiltSpray function. Ignoring the ' ...
%              'FiltSpray function for PIVSurf image']);
%     PIVSurf_Surface_Raw = PIVSurf_Surface_Raw;
%     BadFramePIVSurf=1;
% end
       
 


       
        
       PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Raw, 'robust');
        % Interpolates NaN in the surface
                %plot(PIVSurf_Surface_Int, '.')
        [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 10000); %100000
        %Smoothing of the surface
        if (length(SP.coefs)>length(PIVSurf_Surface_Int))
            PIVSurf_Surface = SP.coefs(2:end-1);
        else
            PIVSurf_Surface=PIVSurf_Surface_Int;
        end
        %PIVSurf_Surface=PIVSurf_Surface; %Correction including vertical offset 
        PIVFused_Surface = PIVSurf_Surface - Y1; 
        % Portion of surface that corresponds to PIV 
        XPIV_PIVSurf_Surface=[X1:1:X1+length(PIVSurf_Surface)-1];
       
        
if(max(PIVFused_Surface) > size(fused_image_rotate_crop,1)+Y1 || min(PIVFused_Surface) < 1)
    BadFramePIVSurf=1;
end