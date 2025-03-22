function [CompVelWater_PIVAir] = transform_PIVSurf_Water_to_PIVAir(T3,Uinv,Vinv,CompVelWater_PIVSurfW,PIV_LFV_Surface,PIV1_A)

%%% Transform from PIVSurf Water to PIV Water
[Xsurf2,Ysurf2] = tformfwd(T3,Uinv,Vinv);
Xsurf = round(Xsurf2(1):Xsurf2(end));
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf);

FIELDS = fieldnames(CompVelWater_PIVSurfW);
for ii =1:length(FIELDS)
    eval(['FIELD = CompVelWater_PIVSurfW.' FIELDS{ii} ';'])
    %% This is the actual resized image; we can retrieve it for checking
    [CompVelWater_PIVA2,Xpos,Ypos] = imtransform(FIELD,T3,'XYScale',1);
    CompVelWater_PIVA2(:,[1:1220,5732:end]) = NaN;
    CompVelWater_PIVA2([1:2560,5889:end],:) = NaN;
    
    % Water velocity with pixel resolution in PIVSurf Water coordinates
    [XX2,YY2] = meshgrid(Xpos(1):Xpos(end),Ypos(1):Ypos(end));
    [XXq2,YYq2] = meshgrid(1:Xpos(end)-1,1:Ypos(end));
    
    CompVelWater_PIVA = interp2(XX2,YY2,CompVelWater_PIVA2,XXq2,YYq2);
    
    %% Warp CompVelWater to match PIV_Surface (PIVSurf Air in PIV Air coordinates)
    movingPoints2 = Ysurf(1103:5710);
    fixedPoints = PIV_LFV_Surface(4333:8940);
    
    %%% Further rototranslation to match PIV_PIVW_Surface and PIV_LFV_Surface
    RotAngle = -12/4608;
    DY2 = 1.5;
    M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
    YY = M*[1:4608;movingPoints2];
    X2 = YY(1,:);
    Y2 = YY(2,:);
    movingPoints = interp1(X2,Y2,1:4608,'linear','extrap');
    
    DeltaY = round(fixedPoints-movingPoints);
    Xmov = 1:size(CompVelWater_PIVA,1);
    CompVelWater_PIVA_warp = nan(size(CompVelWater_PIVA,1),size(CompVelWater_PIVA,2));
    for i = 1:length(DeltaY)
        Xintrp = 1:size(CompVelWater_PIVA,1)/(size(CompVelWater_PIVA,1)-DeltaY(i)):size(CompVelWater_PIVA,1);
        Yintrp = interp1(Xmov,CompVelWater_PIVA(:,i),Xintrp);
        if DeltaY(i) == 0 
            CompVelWater_PIVA_warp(:,i) = CompVelWater_PIVA(:,i);
        elseif DeltaY(i)<0
            CompVelWater_PIVA_warp(:,i) = Yintrp(-DeltaY(i):end);
        else
            CompVelWater_PIVA_warp(:,i) = [nan(1,DeltaY(i)) Yintrp];
        end
    end
    
    %%% Use only masked values
    [PIV_Mask_W1] = PIVWater_Mask(CompVelWater_PIVA, PIV_LFV_Surface(4333:8940));
    CompVelWater_PIVA_warp = CompVelWater_PIVA_warp.*PIV_Mask_W1;
%     CompVelWater_PIVA_warp = CompVelWater_PIVA.*PIV_Mask_W1; % use this only if warping doesn't work
    
    %% CompVelWater in PIVAir coordinates
    eval(['CompVelWater_PIVAir.' FIELDS{ii} ' = CompVelWater_PIVA_warp(4:4:end-4,4:4:end-4);'])
end