function [CompVelWater_PIVAir] = transform_PIVSurf_Water_to_PIVAir(T3,Uinv,Vinv,CompVelWater_PIVSurfW,XPIV_LFV_Surface,PIV_LFV_Surface,PIV1_A)

%%% Transform Surface from PIVSurf Water to PIV Air
[Xsurf2,Ysurf2] = tformfwd(T3,Uinv,Vinv);
Xsurf = round(Xsurf2(1):Xsurf2(end));
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf,'linear','extrap');

FIELDS = fieldnames(CompVelWater_PIVSurfW);
for ii = [1 3] %[1,3,5,8] % Loop only for INTdelx, INTdely, dcor (fitted with zeros)
    eval(['FIELD = CompVelWater_PIVSurfW.' FIELDS{ii} ';'])

    % Fit dcor with zeros
    if strcmp(FIELDS{ii},'dcor')
        FIELD(isnan(FIELD)) = 0;
    end

    %% This is the actual resized image; we can retrieve it for check
    [CompVelWater_PIVA2,Xpos,Ypos] = imtransform(FIELD,T3,'XYScale',1);

    % Water velocity with pixel resolution in PIVSurf Water coordinates
    [XX2,YY2] = meshgrid(Xpos(1):Xpos(end),Ypos(1):Ypos(end));
    [XXq2,YYq2] = meshgrid(-46:4199,1:4628);
    %-46 is the first point of PIVWater in PIVAir x-coordinates that is~NaN
    %4199 is the last point of PIVWater in PIVAir x-coordinates that is~NaN
    %4628 is the first point of PIVWater in PIVAir y-coordinates that is~NaN (from the bottom)

    CompVelWater_PIVA = interp2(XX2,YY2,CompVelWater_PIVA2,XXq2,YYq2);
    CompVelWater_xPIV = XXq2(1,:); % x values of CompVelWater_PIVA in PIVAir coordinates
    CompVelWater_yPIV = YYq2(1,:); % y values of CompVelWater_PIVA in PIVAir coordinates

    %% Warp CompVelWater to match PIV_Surface (PIVSurf Air in PIV Air coordinates)
    % Find common points between PIV_PIVW_Surface (Ysurf) and PIV_LFV_Surface
    I1 = find(XPIV_LFV_Surface == Xsurf(1));
    Iend = find(XPIV_LFV_Surface == Xsurf(end));
    fixedPoints = PIV_LFV_Surface(I1:Iend); %Points that are considered correct to warp CompVelWater

    %%% Further rototranslation to match PIV_PIVW_Surface and PIV_LFV_Surface
    RotAngle2 = (9.12-0.32)/(1-5886);
    DY2 = 9.12;
    M = [cos(RotAngle2) sin(RotAngle2); -sin(RotAngle2) cos(RotAngle2)];
    YY = M*[1:length(Xsurf);Ysurf];
    X2 = YY(1,:);
    Y2 = YY(2,:)-DY2;
    movingPoints = interp1(X2,Y2,1:length(Xsurf),'linear','extrap'); %Points that are moved to warp CompVelWater

    %% Match
    % Find starting points of PIVWater in PIVAir coordinates for surface
    % cropping
    Ix1 = find(Xsurf==CompVelWater_xPIV(1));
    IxEnd = find(Xsurf==CompVelWater_xPIV(end));
    fixedPoints = fixedPoints(Ix1:IxEnd);
    movingPoints = movingPoints(Ix1:IxEnd);

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
    [PIV_Mask_W1] = PIVWater_Mask(CompVelWater_PIVA, fixedPoints);
    CompVelWater_PIVA_warp = CompVelWater_PIVA_warp.*PIV_Mask_W1;
    %     CompVelWater_PIVA_warp = CompVelWater_PIVA.*PIV_Mask_W1; % use this only if warping doesn't work

    %% CompVelWater in PIVAir coordinates
    eval(['CompVelWater_PIVAir.' FIELDS{ii} ' = CompVelWater_PIVA_warp(4:4:end-4,4:4:end-4);'])

    %% Mask not interpolated matrices with NOTnan
% % %     NOTnan = CompVelWater_PIVSurfW.NOTnan;
% % %     [CompVelWater_PIVA2,~,~] = imtransform(NOTnan,T3,'XYScale',1);
% % %     NOTnan = interp2(XX2,YY2,CompVelWater_PIVA2,XXq2,YYq2);
% % %     NOTnan(NOTnan<0.5) = NaN;
% % %     NOTnan(NOTnan>=0.5) = 1;
% % %     if strcmp(FIELDS{ii},'INTdelx')
% % %         eval(['CompVelWater_PIVAir.delta_x = CompVelWater_PIVA_warp(4:4:end-4,4:4:end-4).*NOTnan(4:4:end-4,4:4:end-4);'])
% % %     elseif strcmp(FIELDS{ii},'INTdelz')
% % %         eval(['CompVelWater_PIVAir.delta_z = CompVelWater_PIVA_warp(4:4:end-4,4:4:end-4).*NOTnan(4:4:end-4,4:4:end-4);'])
% % %     elseif strcmp(FIELDS{ii},'dcor')
% % %         eval(['CompVelWater_PIVAir.dcor = CompVelWater_PIVA_warp(4:4:end-4,4:4:end-4).*NOTnan(4:4:end-4,4:4:end-4);'])
% % %         NOTnan(isnan(NOTnan)) = 0;
% % %         eval(['CompVelWater_PIVAir.NOTnan = NOTnan(4:4:end-4,4:4:end-4);'])
% % %     end

end
CompVelWater_PIVAir.xPIV = CompVelWater_xPIV(4:4:end-4);
CompVelWater_PIVAir.yPIV = CompVelWater_yPIV(4:4:end-4);
CompVelWater_PIVAir.Mask = PIV_Mask_W1(4:4:end-4,4:4:end-4);
CompVelWater_PIVAir.PF_surface = fixedPoints(4:4:end-4)/4;
