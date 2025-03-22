function [CompVelWater_PIVSurfW] = transform_PIVWater_to_PIVSurf_Water(CompVelWater,PIV1_W,T2inv)

FIELDS = fieldnames(CompVelWater);
for i = 1:length(FIELDS)
    eval(['FIELD = CompVelWater.' FIELDS{i} ';'])
    % FIELD = CompVelWater.(FIELDS{i});  THIS LINES WORKS ON MATLAB2023a
    % instead of "eval"

    % Loop only for INTdelx, INTdely, dcor (fitted with zeros), Mask
    if i>4 && i<length(FIELDS)
        continue
    elseif strcmp(FIELDS{i},'dcor')
        FIELD(isnan(FIELD)) = 0;
    end
    %% Interpolate PIVWater measurements with pixel resolution
    [XX,YY] = meshgrid(4:4:size(FIELD,2)*4,4:4:size(FIELD,1)*4);
    [XXq,YYq] = meshgrid(1:size(PIV1_W,2),1:size(PIV1_W,1));
    % VVq = griddata(XX,YY,CompVelWater3.INTdelx,XXq,YYq);
    VVq = interp2(XX,YY,FIELD,XXq,YYq);

    %% Rototranslate the image
    % This is a correction due to the rototranslation made after
    % transformation from PIVSurf Water to PIV Water. Note that DY and
    % RotAngle are different due to the different resolution between
    % PIVSurf Water and PIV Water
    DY = 15;
    RotAngle = 4/(3432-806);
    % VVq = VVq.*Mask1_W;
    % VVVq = VVq; VVVq(isnan(VVq)) = 9999;
    IMrot = imrotate(VVq,rad2deg(RotAngle));
    VVVq = IMrot; VVVq(isnan(IMrot)) = 9999;
    IMtr = imtranslate(VVVq,[0 DY]);
    IMtr(IMtr == 0) = NaN;
    IMtr(IMtr==9999) = NaN;

    %% Transform from PIV Water to PIVSurf Water
    [FIELD_inv,XposInv,YposInv] = imtransform(IMtr,T2inv,'XYScale',1);
    FIELD_inv(FIELD_inv == 0) = NaN;

    %% Water velocity with pixel resolution in PIVSurf Water coordinates
    [XX2,YY2] = meshgrid(XposInv(1):XposInv(end),YposInv(1):YposInv(end));
    [XXq2,YYq2] = meshgrid(1:XposInv(end),1:YposInv(end));
    % CompVelWater_PIVSurfW = griddata(XX2,YY2,CompVelWater_inv,XXq2,YYq2);
    FIELD_PIVSurfW = interp2(XX2,YY2,FIELD_inv,XXq2,YYq2);

    eval(['CompVelWater_PIVSurfW.' FIELDS{i} ' = FIELD_PIVSurfW;'])
    % CompVelWater.PIVSurfW.(FIELDS{i}) = FIELD_PIVSurfW;;  THIS LINES WORKS ON MATLAB2023a
    % instead of "eval"

    %% Mask not interpolated matrices with NOTnan
    if ~strcmp(FIELDS{i},'Mask')
        NOTnan = CompVelWater.NOTnan;
        NOTnan = interp2(XX,YY,NOTnan,XXq,YYq);
        IMrot = imrotate(NOTnan,rad2deg(RotAngle));
        VVVq = IMrot; VVVq(isnan(IMrot)) = 9999;
        IMtr = imtranslate(VVVq,[0 DY]);
        IMtr(IMtr==9999) = NaN;
        [NOTnan,XposInv,YposInv] = imtransform(IMtr,T2inv,'XYScale',1);
        [XX2,YY2] = meshgrid(XposInv(1):XposInv(end),YposInv(1):YposInv(end));
        [XXq2,YYq2] = meshgrid(1:XposInv(end),1:YposInv(end));
        NOTnan = interp2(XX2,YY2,NOTnan,XXq2,YYq2);
        NOTnan(NOTnan<0.5) = NaN;
        NOTnan(NOTnan>=0.5) = 1;
        if strcmp(FIELDS{i},'INTdelx')
            eval(['CompVelWater_PIVSurfW.delta_x = FIELD_PIVSurfW.*NOTnan;'])
        elseif strcmp(FIELDS{i},'INTdelz')
            eval(['CompVelWater_PIVSurfW.delta_z = FIELD_PIVSurfW.*NOTnan;'])
        elseif strcmp(FIELDS{i},'dcor')
            eval(['CompVelWater_PIVSurfW.dcor = FIELD_PIVSurfW.*NOTnan;'])
            NOTnan(isnan(NOTnan)) = 0;
            eval(['CompVelWater_PIVSurfW.NOTnan = NOTnan;'])
        end
    end

end