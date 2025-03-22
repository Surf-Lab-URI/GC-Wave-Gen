%%% Test surface for 8 cm depth (tests 06082022)
clear
close all
clc

for ii = 6 % 4:12
    
    Movie = num2str(ii);
    Scene = '2';
    xr = [];
    DM = [];
    M = [];
    DeltaS = [];
    Surf_PIV = [];
    Surf_PIVSURF = [];
    %% Load PIV and reconstruct rectified image
    for i = 0%:2:114
        Str = num2str(i);
        Str2 = num2str(i/2);
        if length(Str)<3
            if length(Str)<2
                Str = ['00' Str];
            else
                Str = ['0' Str];
            end
        end
        if length(Str2)<2
            Str2 = ['0' Str2];
        end

        % For real tests
        Str = ['0' Str];
        Str2 = ['00' Str2];
        try
%             PIV = fliplr(load_Image_IOCoreView_48MP(['\\spray1\d\Shoaling_waves\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\PIV\Movie' Movie '_Scene' Scene '_PIV_' Str '.raw'])); %128 used for filter
            PIV = fliplr(load_Image_IOCoreView_48MP(['E:\Fabio\work\Udel\Fabio\Drag_shallow_water\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\PIV\Movie' Movie '_Scene' Scene '_PIV_' Str '.raw'])); %128 used for filter
        catch
            continue
        end
% % % %         PIV(PIV>300) = 300;
        PIV_lens = CorrectPIVLensDistortion(PIV);
        PIV_CamAngle_Corrected = Correct_Angle_PIV(PIV_lens);
        
% % % %         if ii<13
% % % %             [imSurf_PIV] = FindSurface(PIV_CamAngle_Corrected(5001:6000,:), 5,1);
% % % %             % imSurf_PIV.surface(5468:5730) = NaN;
% % % %             % imSurf_PIV.surface = interp1(1:length(imSurf_PIV.surface),imSurf_PIV.surface,1:length(imSurf_PIV.surface),'pchip','extrap');
% % % %             PIV_Surface_Raw = despike_fab(imSurf_PIV.surface+5000);
% % % %         elseif ii>12 && ii<19
% % % %             [imSurf_PIV] = FindSurface(PIV_CamAngle_Corrected(5001:7000,:), 5,1);
% % % %             % imSurf_PIV.surface(5468:5730) = NaN;
% % % %             % imSurf_PIV.surface = interp1(1:length(imSurf_PIV.surface),imSurf_PIV.surface,1:length(imSurf_PIV.surface),'pchip','extrap');
% % % %             PIV_Surface_Raw = despike_fab(imSurf_PIV.surface+5000);
% % % %         elseif ii>18 && ii<25
% % % %             [imSurf_PIV] = FindSurface(PIV_CamAngle_Corrected(5001:8000,:), 5,1);
% % % %             % imSurf_PIV.surface(5468:5730) = NaN;
% % % %             % imSurf_PIV.surface = interp1(1:length(imSurf_PIV.surface),imSurf_PIV.surface,1:length(imSurf_PIV.surface),'pchip','extrap');
% % % %             PIV_Surface_Raw = despike_fab(imSurf_PIV.surface+5000);
% % % %         end
% % % %             
% % % %         PIV_Surface_Raw = despike_fab(PIV_Surface_Raw);
% % % %         PIV_Surface_Raw = despike_fab(PIV_Surface_Raw);
% % % %         PIV_Surface_Int = filt_spray(PIV_Surface_Raw);
% % % %         PIV_Surface_Int = smoothn(PIV_Surface_Int, 'robust');
% % % %         [SP,~] = spaps(1:length(PIV_Surface_Int), PIV_Surface_Int, 10000); %100000
% % % %         % PIV_Surface = SP.coefs; %%% use for flat surface
% % % %         PIV_Surface = SP.coefs(2:end-1);
% % % %         % PIV_Surface = polyval(polyfit([1,length(imSurf_PIV.surface)],PIV_Surface,1),1:length(imSurf_PIV.surface)); %%% use for flat surface
        
        %% PIVSURF 
%         PIVSURF = fliplr(load_Image_IOCoreView_12MP(['\\spray1\d\Shoaling_waves\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\PIVSURF\Movie' Movie '_Scene' Scene '_PIVSURF_' Str2 '.raw'])); %064 used for filter
        PIVSURF = fliplr(load_Image_IOCoreView_12MP(['E:\Fabio\work\Udel\Fabio\Drag_shallow_water\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\PIVSURF\Movie' Movie '_Scene' Scene '_PIVSURF_' Str2 '.raw'])); %064 used for filter
        %PIVSURF(PIVSURF>100) = 100;
        PIVSURF_lens = CorrectPIVSURFLensDistortion(PIVSURF);
        %%% This polynomial transformation also corrects for residual errors from camera
        %%% distortion
        load PIV_match.mat U1 X1 Xt
        PIV2 = PIVSURF_lens;
        T1 = fitgeotrans(U1,X1,'polynomial',4);
        PIVSURF_CamAngle = PIVSURF_lens;
        PIVSURF_CamAngle(isnan(PIVSURF_CamAngle)) = 0;
        PIVSURF_CamAngle = imwarp(PIVSURF_CamAngle,T1,'linear');
        PIVSURF_CamAngle = imtranslate(PIVSURF_CamAngle, [(X1(1,1)-Xt(1,1)) (X1(1,2)-Xt(1,2))] );
        PIVSURF_CamAngle = fliplr(PIVSURF_CamAngle);
        % % % PIV2(isnan(PIV2)) = 0;
        % % % PIV2_CamAngle_Corrected = PIV2;
        
        %%% Resize PIVSURF on PIV image
        P = load('PIV_PIVsurf_matching_points.mat','Xp','Up','PIV_iP','PIVSURF_iP','PIVSURF_iP2');
        Up2 = [-X1(:,1)+size(PIVSURF_CamAngle,2)+1,X1(:,2)];
        Up2([1:26,end-8*26+1:end],:) = [];
        Up2 = flipud(reshape(permute(reshape(Up2,[26,31,2]),[2,1,3]),26*31,2));
        Up = [ Up2(end-30,:) ; Up2(end,:) ; Up2(31,:) ; Up2(1,:) ]; % The coordinate
        % of four points in the PIV surface image (the grid calibration image).
        Xp = P.Xp; % The coordinates of  the four
        % points in the fused PIV image. The physical location  of these points are
        % the same as PIV surface images. The calibration grid was used to find out
        % the exact same locations for these two images.
        T2_bis = maketform('projective',Up,Xp);
        [Resized_PIVSurf_bis,XPos,YPos] =  imtransform(PIVSURF_CamAngle,T2_bis,'XYScale',1);
        PIVSurf_PIV = Resized_PIVSurf_bis(1:size(PIV_CamAngle_Corrected,1)-YPos(1),1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1));
        
        %%% Find surface on resized PIVSURF
        X = 1:size(PIVSURF_CamAngle,2);
        Xgood = 301:length(X);
        X = Xgood;
        [imSurf_PIVSURF] = FindSurface(PIVSURF_CamAngle,5,1);
        PIVSurf_Surface_Raw = despike_fab(imSurf_PIVSURF.surface(X));
        PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
        PIVSurf_Surface_Raw = despike_fab(PIVSurf_Surface_Raw);
        PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
        PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
        [SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 1000); %100000
        % PIVSurf_Surface = SP.coefs; %%% use for flat surface
        PIVSurf_Surface = SP.coefs(2:end-1);
        % PIVSurf_Surface = polyval(polyfit([X(1),X(end)],PIVSurf_Surface,1),X);  %%% use for flat surface
        PIVSurf_Surface_resized2 = tformfwd([X',PIVSurf_Surface'],T2_bis);
%         PIVSurf_Surface_resized = interp1(PIVSurf_Surface_resized2(:,1),PIVSurf_Surface_resized2(:,2),1:length(PIV_Surface));
        X_PIVSURF2 = PIVSurf_Surface_resized2(:,1);
        PIVSurf_Surface_resized = PIVSurf_Surface_resized2(:,2);
        PIVSurf_SurfMatch_PIV = interp1(PIVSurf_Surface_resized2(:,1),PIVSurf_Surface_resized2(:,2),1:size(PIV_CamAngle_Corrected,2));
        
        %figure;imagesc(Resized_PIVSURF(1:size(PIV_CamAngle_Corrected,1)-YPos(1),1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)));hold on;plot(imSurf_PIVSURF.surface(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)),'r')
        
        %%% Further rotation to correct flat surface
        xrr = PIVSurf_Surface_resized';
        XX = [X_PIVSURF2';xrr];
        Alpha = rad2deg(9.15/5924);
        DX = 25.3894;%34.5595;
        Yref = 5626.3;
        MM = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
        Xr = MM*XX;
%         xr{i+1} = interp1(Xr(1,:),Xr(2,:),1:5924,'pchip','extrap');
        X_PIVSURF = round(X_PIVSURF2(1)):round(X_PIVSURF2(end));
        xr{i+1} = interp1(Xr(1,:),Xr(2,:),X_PIVSURF,'pchip','extrap');
        IndX = [find(X_PIVSURF==1),find(X_PIVSURF==size(PIV_CamAngle_Corrected,2))];
        
        %% LFV
%         LFV = fliplr(load_Image_IOCoreView_12MP(['\\spray1\d\Shoaling_waves\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\LFV\Movie' Movie '_Scene' Scene '_LFV_' Str2 '.raw'])); %064 used for filter
        LFV = fliplr(load_Image_IOCoreView_12MP(['E:\Fabio\work\Udel\Fabio\Drag_shallow_water\data\Tests\Movie' Movie '\Movie' Movie '_Scene' Scene '\RAW\LFV\Movie' Movie '_Scene' Scene '_LFV_' Str2 '.raw'])); %064 used for filter
        LFV_lens = CorrectLFVLensDistortion(LFV);
        LFV_CamAngle = LFV_lens;
        LFV_CamAngle(isnan(LFV_CamAngle)) = 0;
        UU1 = [ 184.5  1457 ; 3788.5 1497 ; 3991.5 100 ; 1.75 100];
        XX1 = [ (184.5+1.75)/2 1477 ; (3788.5+3991.5)/2 1477 ; (3788.5+3991.5)/2 100 ; (184.5+1.75)/2 100 ];
        TT1 = fitgeotrans(UU1,XX1,'projective');
        PIV2_CamAngle_Corrected = imwarp(LFV_CamAngle,TT1,'linear');
        Up22 = [ 965 1026 ; 2947 876.5 ; 3459 1419 ; 968 1417 ];
        Xp22 = [ 1196 1887 ; 7951 1351 ; 9759 3182 ; 1190 3187 ];
        T2_bis2 = maketform('projective',Up22,Xp22);
        [Resized_LFV_bis,XPos2,YPos2] =  imtransform(PIV2_CamAngle_Corrected,T2_bis2,'XYScale',1);
        LFV_PIVSurf = Resized_LFV_bis(1-YPos2(1):7242-YPos2(1),1-XPos2(1):9937-XPos2(1)); % LFV image matching PIVSurf
        
        %%% Find surface on resized PIVSURF
        Y = (1:size(PIV2_CamAngle_Corrected,2));
        Ygood = 351:4170;
        Y = Ygood;
        [imSurf_LFV] = FindSurface(PIV2_CamAngle_Corrected,5,1);
        LFV_Surface_Raw = despike_fab(imSurf_LFV.surface(Y));
        LFV_Surface_Raw = despike_fab(LFV_Surface_Raw);
        LFV_Surface_Raw = despike_fab(LFV_Surface_Raw);
        LFV_Surface_Int = filt_spray(LFV_Surface_Raw);
        LFV_Surface_Int = smoothn(LFV_Surface_Int, 'robust');
        [SP,~] = spaps(1:length(LFV_Surface_Int), LFV_Surface_Int, 1000); %100000
        % PIVSurf_Surface = SP.coefs; %%% use for flat surface
        LFV_Surface = SP.coefs(2:end-1);
        % PIVSurf_Surface = polyval(polyfit([X(1),X(end)],PIVSurf_Surface,1),X);  %%% use for flat surface
        LFV_Surface_resized_PIVSURF2 = tformfwd([Y',LFV_Surface'],T2_bis2);
        X_LFV2 = LFV_Surface_resized_PIVSURF2(:,1)+XPos(1);
        LFV_Surface_resized2 = LFV_Surface_resized_PIVSURF2(:,2)+YPos(1);
        LFV_Surface_resized_PIVSURF = interp1(LFV_Surface_resized_PIVSURF2(:,1),LFV_Surface_resized_PIVSURF2(:,2),1:9937);
        LFV_SurfMatch_PIV = LFV_Surface_resized_PIVSURF(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1))+YPos(1);
        X_LFV = round(X_LFV2(1)):round(X_LFV2(end));
        LFV_Surface_resized = interp1(X_LFV2,LFV_Surface_resized2,X_LFV);
        IndX2 = [find(X_LFV==1),find(X_LFV==size(PIV_CamAngle_Corrected,2))];
        
        %%% Check surface on LFV resized image
% % %         figure;imagesc(Resized_LFV_bis);hold on;plot(LFV_Surface_resized-YPos2(1)-YPos(1),'g')
        
        %%% Check LFV and PIVSURF matching PIV
% % %         figure
% % %         imagesc(LFV_PIVSurf(1:size(PIV_CamAngle_Corrected,1)-YPos(1),1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)))
% % %         hold on
% % %         plot(LFV_Surface_resized_PIVSURF(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)),'g')
% % %         figure
% % %         imagesc(Resized_PIVSurf_bis(1:size(PIV_CamAngle_Corrected,1)-YPos(1),1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)))
% % %         hold on
% % %         plot(PIVSurf_Surface_resized-YPos(1),'b');plot(LFV_Surface_resized_PIVSURF(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)),'g')
% % %         figure
% % %         imagesc(PIV_CamAngle_Corrected(YPos(1):end,:))
% % %         hold on
% % %         plot(PIV_Surface-YPos(1),'r');plot(PIVSurf_Surface_resized-YPos(1),'b');plot(LFV_Surface_resized_PIVSURF(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)),'g')
        
        %%% Further rotation to correct flat surface
        xrr = LFV_Surface_resized;
        XX = [(X_LFV);xrr];
        Alpha = rad2deg(9.15/5924);
        DX2 = 30.3151;
        MM = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
        Xr = MM*XX;
        Xr(:,isnan(Xr(1,:))) = [];
%         xr{i+1} = interp1(Xr(1,:),Xr(2,:),1:5924,'pchip','extrap');
        xr_LFV{i+1} = interp1(Xr(1,:),Xr(2,:),X_LFV,'pchip','extrap');
        % xr2{i+1} = interp1(1:5950/5924:5950,xr{i+1},1:5924);
        figure;imagesc(PIV_CamAngle_Corrected);hold on;
% % % %         plot(PIV_Surface,'r');
        plot(X_PIVSURF,xr{i+1}-DX,'b',X_LFV,xr_LFV{i+1}-DX-DX2,'g')

% % % % % % % % % % % % % % % % % % % % % % % % % %  
        XR = xr{i+1};
        XR_LFV = xr_LFV{i+1};
        eval(['save pair' num2str(i/2) ' XR XR_LFV DX DX2 PIV_CamAngle_Corrected X_PIVSURF X_LFV'])
% % % %             'PIV_Surface ...
% % % %         PIV_Surf{i+1} = PIV_Surface;
        PIVSURF_X{i+1} = X_PIVSURF;
        LFV_X{i+1} = X_LFV;
        PIV_im{i+1} = PIV_CamAngle_Corrected;
        continue
% % % % % % % % % % % % % % % % % % % % % % % % % %

        %% Compare PIVSURF and PIV
        % figure;imagesc(Resized_PIVSURF(1:size(PIV_CamAngle_Corrected,1)-YPos(1),1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)));
        % figure;imagesc(PIV_CamAngle_Corrected(YPos(1):7799,:));hold on;plot(imSurf_PIVSURF.surface(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1)),'r')
        M(i+1) = mean(PIVSurf_Surface_resized-PIV_Surface);
        DeltaS{i+1} = PIVSurf_Surface_resized-PIV_Surface-M(i+1);
        % % %     figure;plot(PIVSurf_Surface_resized,'r');hold on;plot(PIV_Surface,'b')
        % % %     figure;imagesc(PIV_CamAngle_Corrected);hold on;plot(PIVSurf_Surface_resized,'r');plot(PIV_Surface,'b')
        % % %     figure;plot(DeltaS{i},'b')
        
        Surf_PIV{i+1} = PIV_Surface;
        Surf_PIVSURF{i+1} = PIVSurf_Surface_resized;
        PIV_im{i+1} = PIV_CamAngle_Corrected;
        PIVSURF_im{i+1} = Resized_PIVSurf_bis;
        DM{i+1} = mean(xr{i+1}-Surf_PIV{i+1});
        
        %% Reconstruct final surface from PIVSURF to PIV
        %%% Use this section to build a filter for each image %%%
        % x = imSurf_PIVSURF.surface(1-XPos(1):size(PIV_CamAngle_Corrected,2)-XPos(1));
        % y = imSurf_PIV.surface-YPos(1)+1;
        % figure; plot(x-mean(x-y,'omitnan'),'r');hold on;plot(y,'b')
        % [Surface] = Reconstruct_Surface_h8cm(x,y);
        % figure;imagesc(PIV_CamAngle_Corrected)
        % hold on
        % plot(imSurf_PIV.surface,'r')
        % plot(Surface+YPos(1),'k')
        
        %% Use filter based on 05272022_Movie1_Scene2_PIV_128 (064 PIVSURF)
        % x = PIVSurf_Surface_resized(1:end-80);
        % y = PIV_Surface(81:end);
        % figure; plot(x-mean(x-y,'omitnan'),'r');hold on;plot(y,'b')
        % [Surface] = Surface_from_PIVSURF_to_PIV_h8cm(x);
        % figure;imagesc(PIV_CamAngle_Corrected)
        % hold on
        % plot(PIV_Surface,'r')
        % plot(81:length(PIV_Surface),Surface,'k')
        
    end
    
% % % %     eval(['Test' num2str(ii) '.xr = xr;'])
% % % %     eval(['Test' num2str(ii) '.DM = DM;'])
% % % %     eval(['Test' num2str(ii) '.M = M;'])
% % % %     eval(['Test' num2str(ii) '.DeltaS = DeltaS;'])
% % % %     eval(['Test' num2str(ii) '.Surf_PIV = Surf_PIV;'])
% % % %     eval(['Test' num2str(ii) '.Surf_PIVSURF = Surf_PIVSURF;'])
% % % %     
% % % %     eval(['save check_surf Test' num2str(ii) ' -append'])
% % % %     eval(['clear Test' num2str(ii) ])
end