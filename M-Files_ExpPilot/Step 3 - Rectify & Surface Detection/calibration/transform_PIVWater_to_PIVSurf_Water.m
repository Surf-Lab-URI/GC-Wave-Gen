function [CompVelWater_PIVSurfW] = transform_PIVWater_to_PIVSurf_Water(CompVelWater,PIV1_W,T2inv)

FIELDS = fieldnames(CompVelWater);
for i = 1:length(FIELDS)
    eval(['FIELD = CompVelWater.' FIELDS{i} ';'])
    if size(FIELD,1)<2 || size(FIELD,2)<2
        continue
    end
    %% Interpolate PIVWater measurements with pixel resolution
    [XX,YY] = meshgrid(4:4:size(FIELD,2)*4,4:4:size(FIELD,1)*4);
    [XXq,YYq] = meshgrid(1:size(PIV1_W,2),1:size(PIV1_W,1));
    % VVq = griddata(XX,YY,CompVelWater3.INTdelx,XXq,YYq);
    VVq = interp2(XX,YY,FIELD,XXq,YYq);
    
    %% Rototranslate the image
    DY = 17; % Same translation used in Extract_PIVSurf Water
    RotAngle = -39/6107; % Same rotation angle used in Extract_PIVSurf Water
    Minv = [cos(-RotAngle) sin(-RotAngle); -sin(-RotAngle) cos(-RotAngle)];
    
    VVq = VVq; %.*Mask1_W;
    % VVVq = VVq; VVVq(isnan(VVq)) = 9999;
    IMrot = imrotate(VVq,-rad2deg(RotAngle));
    VVVq = IMrot; VVVq(isnan(IMrot)) = 9999;
    IMtr = imtranslate(VVVq,[0 -45]);
    IMtr(IMtr==9999) = NaN;
    %IMrot = imrotate(IMtr,-rad2deg(RotAngle));
    
    %% Transform from PIV Water to PIVSurf Water
    % [CompVelW_inv,Xpos,Ypos] = imtransform(IMrot,T2inv,'XYScale',1);
    [FIELD_inv,XposInv,YposInv] = imtransform(IMtr,T2inv,'XYScale',1);
    
    %% Water velocity with pixel resolution in PIVSurf Water coordinates
    [XX2,YY2] = meshgrid(XposInv(1):XposInv(end),YposInv(1):YposInv(end));
    [XXq2,YYq2] = meshgrid(1:XposInv(end),1:YposInv(end));
    % CompVelWater_PIVSurfW = griddata(XX2,YY2,CompVelWater_inv,XXq2,YYq2);
    FIELD_PIVSurfW = interp2(XX2,YY2,FIELD_inv,XXq2,YYq2);
    eval(['CompVelWater_PIVSurfW.' FIELDS{i} ' = FIELD_PIVSurfW;'])
end