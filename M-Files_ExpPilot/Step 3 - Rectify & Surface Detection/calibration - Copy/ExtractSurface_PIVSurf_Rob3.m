function [BadFramePIVSurf,XPIV_PIVSurf_Surface,PIVSurf_Surface,PIVFused_Surface,PIVSurf_PIVMatch, index] = ExtractSurface_PIVSurf_Rob3(norm_ptv_surf,PIVSurf_CRR,fused_image_rotate_crop,ptv_surf_guide )

avg = nanmean(nanmedian(PIVSurf_CRR.img(end-3000:end-300, 200:end-200)));
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
% PIVSurf_PIVMatch_surf_test = zeros(size(PIVSurf_PIVMatch));
% 
% for i = 1:length(PIVSurf_PIVMatch(1,:))
%  index5 = find(abs(diff(PIVSurf_PIVMatch(:, i)))>0.5);
%  PIVSurf_PIVMatch_surf_test(index5+1, i) = 1;
% 
% end
% figure()
% imagesc(PIVSurf_PIVMatch_surf_test)
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

% PIVSurf image exactely matching height of PIV images, but longer
% end
% index = nanmin(find(Surf == nanmax(Surf)))

% NN=360;
% PIVSurf_Strip = PIVSurf_CRR.img(Y1+1500:Y2, 1+NN:end-NN);  %Correct size for matlabe 2015 and later
% PIVSurf_Strip=PIVSurf_Strip./NORM_PIV_SURF;
% PIVSurf_Strip=imfilter(PIVSurf_Strip,h,'replicate');
% 361 gets rid of NaN from the rotation and matching of PIVSurf with PIV images
%1500 is a vertical offset to get closer to the surface and reduce image
%size in which to detect the surface

%PF_XPos_in_PIVSurf = X1-360:X2-360; %to re-place surface from strip into PIV coordinates

%Creating filter around surface
%PIVSurf_PIVMatch=imgaussian(PIVSurf_PIVMatch,200);
% index = nanmean(ptv_surf_guide);
index = 7500;
PIVSurf_Strip = PIVSurf_PIVMatch(nanmin(nanmax(index-4100, 101),length(PIVSurf_PIVMatch(:,1)) -8100) :nanmin(nanmax(index-4100, 101)+8000,length(PIVSurf_PIVMatch(:,1))-100) , :);
PIVSurf_Strip_T = PIVSurf_PIVMatch(nanmin(nanmax(index-4000, 201),length(PIVSurf_PIVMatch(:,1)) -8000) :nanmin(nanmax(index-4000, 201)+8000,length(PIVSurf_PIVMatch(:,1))) , :);
PIVSurf_Strip_B = PIVSurf_PIVMatch(nanmin(nanmax(index-4200, 1),length(PIVSurf_PIVMatch(:,1)) -8200) :nanmin(nanmax(index-4200, 1)+8000,length(PIVSurf_PIVMatch(:,1))-200) , :);

%PTV_SURF_GUIDE_CORRECTED = ptv_surf_guide - nanmin(nanmax(index-4100, 101),length(PIVSurf_PIVMatch(:,1) -8100)) -100;
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
%hold on 
%plot(PTV_SURF_GUIDE_CORRECTED, '.')
S_positive = S;
S_negative = S;
S_positive(S_positive<0) = 0;
S_negative(S_negative>0) = 0;

% imshow(S_negative, [])
% imshow(S_positive, [])
S_check = [zeros(200,3931); S_negative(200:end-1, :) + S_positive(1:end-200, :)];
imshow(S_check, [])
[BW5, L, N, A] = bwboundaries(imbinarize(S_check, 8));

%%
bw_6 = zeros(size(S_check));.
.
.
.

for i = 1:length(BW5)
    
     plot(BW5{i,1}(:,2),BW5{i,1}(:,1), '.')
     hold on 
  
end

index10 = find(cellfun('length', BW5)>1000)
for f = 1:length(BW5)
    
  if (ismember(f,index10))
    for L =nanmin(BW5{f,1}(:,2)):nanmax(BW5{f,1}(:,2))
        L;
      indextest = [];
      indextest =  find(BW5{f,1}(:,2) == L);
      try
          range(BW5{f,1}(indextest,1))
      if (range(BW5{f,1}(indextest,1)) <300)
          if (range(BW5{f,1}(indextest,1)) > 100)
             bw_6(nanmax(BW5{f,1}(indextest,1)),L) = 1;
             HEIGHT42 = nanmax(BW5{f,1}(indextest,1));
          else
          end
      end
      catch
      end
    end 
  else
  end
end
  imshow(bw_6, [])
  SE = strel('rectangle',[50, 5])
  bw_6thick = imdilate(bw_6, se);




 
  PIVSurf_Surface_Tem = FindSurface(S_check,1,bw_6thick);
  PIVSurf_Surface_Tem.surface =  PIVSurf_Surface_Tem.surface;
imshow(S_check, [])
hold on 
plot(PIVSurf_Surface_Tem.surface, '.r')
f
%%


% index2 = find(cellfun('length', BW5) == nanmax(cellfun('length', BW5)));
% 
% bw_6 = zeros(size(S_check));
% figure
% for i = 1:length(index2)
%      plot(BW5{i,1}(:,2),BW5{i,1}(:,1), '.')
%      hold on 
%     for h = 1:length(BW5{index2(i),1})
%          bw_6(BW5{index2(i),1}(h,1),BW5{index2(i),1}(h,2)) = 1;
%     end
% end
% PIVSurf_Surface_Tem = FindSurface(S_check,1,bw_6);
% imshow(S_check, [])
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.r')
% var_test = zeros(size(S_check));
% for h = nanmax(nanmin(PIVSurf_Surface_Tem_estimate.surface)-300, 2):nanmin(nanmax(PIVSurf_Surface_Tem_estimate.surface)+300, length(S_check(:,1))-1)
%     h
%     for g = 2:length(S_check(1,:))-1
%         var_test(h,g) = nanvar([S_check(h-1:h+1, g-1); S_check(h-1:h+1, g); S_check(h-1:h+1, g+1)]);
%     end
% end
% imshow(var_test, [])
% % [BW5, L, N, A] = bwboundaries(imbinarize(var_test, 0.1));
% % ho
%% imshow(S_check, [])
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.')
% S_check = S_check+[zeros(100,3931); S_negative(99:end-2, :) + S_positive(1:end-100, :)];
% S_check = S_check+[zeros(101,3931); S_negative(10:end-1, :) + S_positive(1:end-101, :)];
% S_check = S_check+[zeros(100,3931); S_negative(100:end-1, :) + S_positive(1:end-100, :)];
% figure(3)
% imshow, [])
% SE = strel('disk',300);
% S4 =  imdilate(S_check,SE);
% imshow(S3, [])
% 
% 
% 
% S4 = edge(S_check, 'prewitt');
%  S3_test = zeros(size(S3));
%  for i = 1:length(S3(1,:))
%     
%     index5 = find((abs(diff((S3(:, i)))))>0.25);
%     S3_test(index5+1, i) = 1;
% 
% end
% figure()
% imagesc(S3_test)
%[BW5, L, N, A] = bwboundaries(imbinarize(S4, 12));
% 

%     %PIVSurf_Strip(BW5{index2,1}(h,1),BW5{d,1}(h,2))
%     for h = 1:length(BW5{d,1})
%         POINTS(h) =  PIVSurf_Strip(BW5{d,1}(h,1),BW5{d,1}(h,2));
%     end
%     BOUNDARY_VAL(d) = nanmean(POINTS);
%     



    
%end
% 
% index2 = find(abs(BOUNDARY_VAL - avg)<20);
% bw_6 = zeros(size(S_check));
% figure
% for i = 1:length(index2)
%      plot(BW5{i,1}(:,2),BW5{i,1}(:,1), '.')
%      hold on 
%     for h = 1:length(BW5{index2(i),1})
%          BW5{index2(i),1}(h,2)
%          BW5{index2(i),1}(h,1)
%          bw_6(BW5{index2(i),1}(h,2),BW5{index2(i),1}(h,1)) = 1;
%     end
% end
% f = fit(X_check', Y_check', 'poly2')
% param.minMajorAxis = 1;
% param.maxMajorAxis = 10000;
% param.rotation = 0;
% param.rotationSpan = 180;
% param.minAspectRatio =  .0001;
% param.numBest = 1; 
% param.randomize = 3; 
%     
% bestFits_guess_surf = ellipseDetection(bw_6,param);
% coeff_fit = coeffvalues(f);
%  for i = 1:3931
%      guess_line(i) =coeff_fit(3) + coeff_fit(2)*i + coeff_fit(1).*(i^(2));
%  end

% for h = 1:length(BW5)
%    if (nanmin(abs(h - index2)) == 0)
%   plot(BW5{h,1}(:,2),BW5{h,1}(:,1), '.')
%   hold on
%    else
%    end
% end 
% plot(guess_line, '.k')
% clear index2
% for g = 1:length(PTV_SURF_GUIDE_CORRECTED)
%     for h = 1:length(BW5)
%        dis_part(h) = nanmin(sqrt((g-BW5{h,1}(:,2)).^(2) + (PTV_SURF_GUIDE_CORRECTED(g) - BW5{h,1}(:,1)).^(2)));
%     end
%     index2(g) = find(dis_part ==nanmin(dis_part));
% end
%%

% 
% index2 = find(cellfun('length', BW5) == nanmax(cellfun('length', BW5)));
% 
% 
% 
% 
% 
% if ((range(BW5{index2, 1}(:,2))./3931) >0.95)
%     D = size(S);
%     S_BW = zeros(D(1), D(2));
%     for L = 1:length(BW5)
%         if (nanmin(abs(L - index2)) == 0)
%             for k = 1:cellfun('length', BW5(L))
%                 S_BW(BW5{L,1}(k,1), BW5{L,1}(k,2)) = 1;
%                 X(k) = BW5{L,1}(k,2);
%                 Y(k) = BW5{L,1}(k,1);
%             end
%             
%         end 
%     end
%     SE = strel('disk',300);
%     S_BW2 =  imdilate(S_BW,SE);
%     S_BW3 = imfill(S_BW2);
%      S3 = (S_check).*S_BW3;
%      imshow(S3, [])
%      PIVSurf_Surface_Tem = FindSurface(PIVSurf_Strip,5,S3); %5
%      plot(PIVSurf_Surface_Tem.surface, '.')
%      index_filt = find(abs(diff(PIVSurf_Surface_Tem.surface)) >20);
% 
%      difference = diff(PIVSurf_Surface_Tem.surface);
%      for f =1:length(difference)
%         inte(f) =nansum(difference(1:f));
%      end
%    
%      %plot(inte, '.')
% 
% 
%      up = find(difference>20);
%      down = find(difference<-20);
%      up = up+1;
%      down = down+1;
% 
%      length(up)
%      length(down)
% 
%      plot(up,PIVSurf_Surface_Tem.surface(up), 'o')
%      hold on 
%      plot(down,PIVSurf_Surface_Tem.surface(down), 'o')
% 
% Test = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Test(down) = 0;
% Test(up) = 1;
% Surface1 = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Surface2 = NaN(1, length(PIVSurf_Surface_Tem.surface));
% 
% 
% index_1 = find(isnan(Test) == 0);
% if (Test(index_1(1)) ==1)
%     Surface1(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1);
%     check = 0
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%             Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
%     
% else
%     Surface2(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1);
%     check = 1 
%    
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
% 
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%              Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
% % end
% end
%    
%      surf1 = length(find(isnan(Surface1)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf1));
%      surf2 = length(find(isnan(Surface2)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf2));
%            if (surf1 < surf2)
%                
%               PIVSurf_Surface_Tem.surface = Surface2;
%               PIVSurf_Surface_Int_surf2 = smoothn(Surface2, 'robust');
%            else
%               PIVSurf_Surface_Tem.surface = Surface1;
%               PIVSurf_Surface_Int_surf1 = smoothn(Surface1, 'robust');
%            end
% else
%      S3 = S4;
%      PIVSurf_Surface_Tem = FindSurface(PIVSurf_Strip,4,S3); %5
%      
%      index_filt = find(abs(diff(PIVSurf_Surface_Tem.surface)) >20);
% 
%      difference = diff(PIVSurf_Surface_Tem.surface);
%      for f =1:length(difference)
%         inte(f) =nansum(difference(1:f));
%      end
%      plot(inte, '.')
%      
% 
%      up = find(difference>20);
%      down = find(difference<-20);
%      up = up+1;
%      down = down+1;
% 
%      length(up)
%      length(down)
% 
%      plot(up,PIVSurf_Surface_Tem.surface(up), 'o')
%      hold on 
%      plot(down,PIVSurf_Surface_Tem.surface(down), 'o')
% 
% Test = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Test(down) = 0;
% Test(up) = 1;
% Surface1 = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Surface2 = NaN(1, length(PIVSurf_Surface_Tem.surface));
% 
% 
% index_1 = find(isnan(Test) == 0);
% if (Test(index_1(1)) ==1)
%     Surface1(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1);
%     check = 0
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%             Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
%     
% else
%     Surface2(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1);
%     check = 1;
%    
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
% 
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%              Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
% % end
% end
%    
%      surf1 = length(find(isnan(Surface1)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf1));
%      surf2 = length(find(isnan(Surface2)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf2));
%            if (surf1 < surf2)
%                
%               PIVSurf_Surface_Tem.surface = Surface2;
%               PIVSurf_Surface_Int_surf2 = smoothn(Surface2, 'robust');
%            else
%               PIVSurf_Surface_Tem.surface = Surface1;
%               PIVSurf_Surface_Int_surf1 = smoothn(Surface1, 'robust');
%            end
% end

%%

% %  %%
%  imshow(S2, [])
%  hold on
%  SE = strel('disk',40);
% S3 =  imdilate(S2,SE);

%PIVSurf_Surface_Tem = FindSurface(PIVSurf_Strip,4,S3); %5
%plot(PIVSurf_Surface_Tem.surface, '.')
% index_filt = find(abs(diff(PIVSurf_Surface_Tem.surface)) >200);
% 
% difference = diff(PIVSurf_Surface_Tem.surface);
% for f =1:length(difference)
%     inte(f) =nansum(difference(1:f));
% end
% plot(inte, '.')
% 
% 
% up = find(difference>200);
% down = find(difference<-200);
% up = up+1;
% down = down+1;
% 
% length(up)
% length(down)
% 
% plot(up,PIVSurf_Surface_Tem.surface(up), 'o')
% hold on 
% plot(down,PIVSurf_Surface_Tem.surface(down), 'o')


% Test = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Test(down) = 0;
% Test(up) = 1;
% Surface1 = NaN(1, length(PIVSurf_Surface_Tem.surface));
% Surface2 = NaN(1, length(PIVSurf_Surface_Tem.surface));


% index_1 = find(isnan(Test) == 0)
% if (Test(index_1(1)) ==1)
%     Surface1(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1)
%     check = 0
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%             Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
%     
% else
%     Surface2(1:index_1(1)-1) = PIVSurf_Surface_Tem.surface(1:index_1(1)-1)
%     check = 1 
%    
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
% 
%         elseif (check ==0)
%              Surface2(j) = PIVSurf_Surface_Tem.surface(j);
%         else
%              Surface1(j) = PIVSurf_Surface_Tem.surface(j);
%         end 
%     end
% % end
%   PIVSurf_Surface_Tem.surface(up) = NaN;
%               PIVSurf_Surface_Tem.surface(down) = NaN;
% %                 Surface2(up) = NaN;
% %                 Surface2(down) = NaN;
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.b')
% smoothsurf = movavg(PIVSurf_Surface_Tem.surface', 'simple', 50);
% hold on 



%  index3 = find(abs(PIVSurf_Surface_Tem.surface' - smoothsurf) > 75)
% PIVSurf_Surface_Tem.surface(index3) = NaN;
% plot(PIVSurf_Surface_Tem.surface, '.r')
%plot(Surface2, '.r')


%    surf1 = length(find(isnan(Surface1)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf1));
%      surf2 = length(find(isnan(Surface2)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf2));
%            if (surf1 < surf2)
%                
%                  PIVSurf_Surface_Tem.surface = Surface2;
% %              PIVSurf_Surface_Int_surf2 = smoothn(Surface2, 'robust');
% %              
% %                
% %                                     figure() 
% %                            imshow(S2, [])
% %                      hold on
% %                           index_filt = find(abs(diff(PIVSurf_Surface_Int_surf2)) >10);
% % 
% %                     difference = diff(PIVSurf_Surface_Int_surf2);
% %                     up = find(difference>10);
% %                     down = find(difference<-10);
% %                     up = up+1;
% %                     down = down+1;
% % 
% %                     length(up)
% %                     length(down)
% % 
% %                     plot(up,PIVSurf_Surface_Int_surf2(up), 'o')
% %                     hold on 
% %                     plot(down,PIVSurf_Surface_Int_surf2(down), 'o')
% % 
% % 
% %                     Test = NaN(1, length(PIVSurf_Surface_Int_surf2));
% %                     Test(down) = 0;
% %                     Test(up) = 1;
% %                     Surface3 = NaN(1, length(PIVSurf_Surface_Int_surf2));
% %                     Surface4 = NaN(1, length(PIVSurf_Surface_Int_surf2));
% % 
% %                     index_1 = find(isnan(Test) == 0)
% %                     if (Test(index_1(1)) ==1)
% %                         Surface3(1:index_1(1)-1) = PIVSurf_Surface_Int_surf2(1:index_1(1)-1)
% %                         check = 0
% %                         for j = index_1(1):length(Test)
% %                             if (Test(j) == check)
% %                                 check = check+1;
% %                                 check = mod(check, 2)
% %                             elseif (check ==0)
% %                                  Surface4(j) =  PIVSurf_Surface_Int_surf2(j);
% %                             else
% %                                 Surface3(j) =  PIVSurf_Surface_Int_surf2(j);
% %                             end 
% %                         end
% % 
% %                     else
% %                         Surface4(1:index_1(1)-1) =  PIVSurf_Surface_Int_surf2(1:index_1(1)-1)
% %                         check = 1 
% % 
% %                         for j = index_1(1):length(Test)
% %                             if (Test(j) == check)
% %                                 check = check+1;
% %                                 check = mod(check, 2)
% % 
% %                             elseif (check ==0)
% %                                  Surface4(j) =  PIVSurf_Surface_Int_surf2(j);
% %                             else
% %                                 Surface3(j) =  PIVSurf_Surface_Int_surf2(j);
% %                             end 
% %                         end
% %                     end
% %                       Surface3(up) = NaN;
% %                 Surface3(down) = NaN;
% %                 Surface4(up) = NaN;
% %                 Surface4(down) = NaN;
% %                     hold on 
% %                     plot(Surface3, '.b')
% %                     hold on 
% %                     plot(Surface4, '.r')
% % 
% %              surf3 = length(find(isnan(Surface3)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf1));
% %              surf4 = length(find(isnan(Surface4)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf2));
% %            if (surf3 < surf4)
% %                PIVSurf_Surface_Tem.surface = Surface4;
% %            elseif (surf3 > surf4)
% %                 PIVSurf_Surface_Tem.surface = Surface3;
% %            else
% %                 PIVSurf_Surface_Tem.surface = (Surface4+Surface3)./2;
% %            end 
% 
%        elseif (surf1 > surf2)
%             PIVSurf_Surface_Tem.surface = Surface1;
% %              PIVSurf_Surface_Int_surf1 = smoothn(Surface1, 'robust');
% %                                 figure() 
% %                        imshow(S2, [])
% %                  hold on
% %                       index_filt = find(abs(diff(PIVSurf_Surface_Int_surf1)) >10);
% % 
% %                 difference = diff(PIVSurf_Surface_Int_surf1);
% %                 up = find(difference>10);
% %                 down = find(difference<-10);
% %                 up = up+1;
% %                 down = down+1;
% % 
% %                 length(up)
% %                 length(down)
% % 
% %                 plot(up,PIVSurf_Surface_Int_surf1(up), 'o')
% %                 hold on 
% %                 plot(down,PIVSurf_Surface_Int_surf1(down), 'o')
% % 
% % 
% %                 Test = NaN(1, length(PIVSurf_Surface_Int_surf1));
% %                 Test(down) = 0;
% %                 Test(up) = 1;
% %                 Surface3 = NaN(1, length(PIVSurf_Surface_Int_surf1));
% %                 Surface4 = NaN(1, length(PIVSurf_Surface_Int_surf1));
% % 
% %                 index_1 = find(isnan(Test) == 0)
% %                 if (Test(index_1(1)) ==1)
% %                     Surface3(1:index_1(1)-1) = PIVSurf_Surface_Int_surf1(1:index_1(1)-1)
% %                     check = 0
% %                     for j = index_1(1):length(Test)
% %                         if (Test(j) == check)
% %                             check = check+1;
% %                             check = mod(check, 2)
% %                         elseif (check ==0)
% %                              Surface4(j) =  PIVSurf_Surface_Int_surf1(j);
% %                         else
% %                             Surface3(j) =  PIVSurf_Surface_Int_surf1(j);
% %                         end 
% %                     end
% % 
% %                 else
% %                     Surface4(1:index_1(1)-1) =  PIVSurf_Surface_Int_surf1(1:index_1(1)-1)
% %                     check = 1 
% % 
% %                     for j = index_1(1):length(Test)
% %                         if (Test(j) == check)
% %                             check = check+1;
% %                             check = mod(check, 2)
% % 
% %                         elseif (check ==0)
% %                              Surface4(j) =  PIVSurf_Surface_Int_surf1(j);
% %                         else
% %                             Surface3(j) =  PIVSurf_Surface_Int_surf1(j);
% %                         end 
% %                     end
% %                 end
% %                 Surface3(up) = NaN;
% %                 Surface3(down) = NaN;
% %                 Surface4(up) = NaN;
% %                 Surface4(down) = NaN;
% %                 hold on 
% %                 plot(Surface3, '.b')
% %                 hold on 
% %                 plot(Surface4, '.r')
% % 
% %          surf3 = length(find(isnan(Surface3)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf1));
% %          surf4 = length(find(isnan(Surface4)==0));         %nansum(abs(PTV_SURF_GUIDE_CORRECTED - PIVSurf_Surface_Int_surf2));
% %            if (surf3 < surf4)
% %                PIVSurf_Surface_Tem.surface = Surface4;
% %            elseif (surf3 > surf4)
% %                 PIVSurf_Surface_Tem.surface = Surface3;
% %            else
% %                 PIVSurf_Surface_Tem.surface = (Surface4+Surface3)./2;
% %            end 
%            
%            
%            
%            
%            
%            else
%                PIVSurf_Surface_Tem.surface = Surface1;
%        end

                %    
 
   %%
          
       
%       index_filt = find(abs(diff(PIVSurf_Surface_Int_surf2)) >20);
% 
% difference = diff(PIVSurf_Surface_Int_surf2);
% up = find(difference>20);
% down = find(difference<-20);
% up = up+1;
% down = down+1;
% 
% length(up)
% length(down)
% 
% plot(up,PIVSurf_Surface_Int_surf2(up), 'o')
% hold on 
% plot(down,PIVSurf_Surface_Int_surf2(down), 'o')
% 
% 
% Test = NaN(1, length(PIVSurf_Surface_Int_surf2));
% Test(down) = 0;
% Test(up) = 1;
% Surface5 = NaN(1, length(PIVSurf_Surface_Int_surf2));
% Surface6 = NaN(1, length(PIVSurf_Surface_Int_surf2));
% 
% index_1 = find(isnan(Test) == 0)
% if (Test(index_1(1)) ==1)
%     Surface5(1:index_1(1)-1) = PIVSurf_Surface_Int_surf2(1:index_1(1)-1)
%     check = 0
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
%         elseif (check ==0)
%              Surface6(j) =  PIVSurf_Surface_Int_surf2(j);
%         else
%             Surface5(j) =  PIVSurf_Surface_Int_surf2(j);
%         end 
%     end
%     
% else
%     Surface6(1:index_1(1)-1) =  PIVSurf_Surface_Int_surf2(1:index_1(1)-1)
%     check = 1 
%    
%     for j = index_1(1):length(Test)
%         if (Test(j) == check)
%             check = check+1;
%             check = mod(check, 2)
% 
%         elseif (check ==0)
%              Surface6(j) =  PIVSurf_Surface_Int_surf2(j);
%         else
%             Surface5(j) =  PIVSurf_Surface_Int_surf2(j);
%         end 
%     end
% end
% hold on 
% plot(Surface5, '.b')
% hold on 
% plot(Surface6, '.r')
%    
   
       
 %%      
       
       


%%
%PIVSurf_Surface_Tem.surface(index_filt+1) = NaN;
%plot(PIVSurf_Surface_Tem.surface, '.')
% [fitresult, gof] = nyqist_guess_rob(find(isnan(PIVSurf_Surface_Tem.surface) == 0), PIVSurf_Surface_Tem.surface(isnan(PIVSurf_Surface_Tem.surface) == 0));
% if (gof.rsquare >0.35)
%      params =  coeffvalues(fitresult);
%      phase_longwave_ptv_surf = params(3);
%      wavenum_longwave = params(2);
%      offset = params(4);
%      wave_height = params(1);
%      SURF_TEST = PIVSurf_Surface_Tem.surface;
%      SURF_TEST(isnan(PIVSurf_Surface_Tem.surface) == 0) = wave_height*cos(wavenum_longwave*PIVSurf_Surface_Tem.surface(isnan(PIVSurf_Surface_Tem.surface) == 0) +  phase_longwave_ptv_surf) + offset;
%      PIVSurf_Surface_Tem.surface = SURF_TEST;
% else
%finds surface in the strip
% imshow(S, [])
% hold on 
% plot(PIVSurf_Surface_Tem.surface, '.')
%index_NaN = (find(isnan(PIVSurf_Surface_Tem.surface))); 
% gradient_check = nanmean(diff(PIVSurf_Surface_Tem.surface));
% gradient2_check = nanmedian(diff(diff(PIVSurf_Surface_Tem.surface))); 
%gradient3_check = nanmean(diff(diff(diff(PIVSurf_Surface_Tem.surface)))); 

% SURF_TEST = PIVSurf_Surface_Tem.surface;
% try
% for j = 1:length(index_NaN)
%    
%     try
%         SURF_TEST(index_NaN(j)) = gradient_check*(index_NaN(j)-(nanmax(index_NaN+1)))+ (gradient2_check./2).*((index_NaN(j)- (nanmax(index_NaN+1)))^2) + SURF_TEST((nanmax(index_NaN+1)));  %+ (gradient2_check./2).*((index_NaN(j)- (nanmax(index_NaN+1)))^2) +  (gradient3_check./6).*((index_NaN(j)- (nanmax(index_NaN+1)))^3)
%     catch
%         SURF_TEST(index_NaN(j)) = gradient_check*(index_NaN(j)-(nanmin(index_NaN-1))) + (gradient2_check./2).*((index_NaN(j)- (nanmin(index_NaN-1)))^2) + SURF_TEST((nanmin(index_NaN-1))); %(gradient2_check./2).*((index_NaN(j)- (nanmin(index_NaN-1)))^2) + (gradient3_check./6).*((index_NaN(j)- (nanmin(index_NaN-1)))^3)+ 
%     end  
% end 
%PIVSurf_Surface_Tem.surface(index_NaN) = PTV_SURF_GUIDE_CORRECTED(index_NaN);
% catch
% end


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