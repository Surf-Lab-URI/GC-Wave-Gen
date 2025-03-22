function [Cartesian] = RemoveOutliers_fabio(CompVel,Threshold)

delx=CompVel.delta_x;
delz=CompVel.delta_z;
dcor=CompVel.dcor;

% Find bad correlation
delx(dcor < Threshold) = NaN;
delz(dcor < Threshold) = NaN;
nc=sum(sum((dcor < Threshold)));

% Size of figure, NaN in Mask, and total NaN
s=size(dcor); s=s(1)*s(2);
s_Mask = sum(isnan(CompVel.Mask),'all');
s_NaN = sum(isnan(dcor),'all');
delx(isnan(dcor)) = NaN;
delz(isnan(dcor)) = NaN;

% Find bad subpixel velocities
ThrX = (delx-CompVel.INTdelx)>3 | (delx-CompVel.INTdelx)<-3;
ThrZ = (delz-CompVel.INTdelz)>2 | (delz-CompVel.INTdelz)<-2;
delx(ThrX) = NaN; % NaN outliers
delz(ThrZ) = NaN; % NaN outliers
np = sum([ sum(abs(delx-CompVel.INTdelx)>3,'all') , sum(abs(delz-CompVel.INTdelz)>2,'all')]);


%% Outlier removed
Cartesian.delx=delx.*CompVel.Mask;%.*MaskNaN; if dcor and delx has same NaN
Cartesian.delz=delz.*CompVel.Mask;%.*MaskNaN; if dcor and delx has same NaN

%% Smoothed
Cartesian.u=smoothn(Cartesian.delx,0.4,'robust');
Cartesian.w=smoothn(Cartesian.delz,0.4,'robust');
Cartesian.corr_outlier_percent=nc/(s-s_Mask)*100;
Cartesian.subpix_outlier_percent=np/(s-s_Mask)*100;
Cartesian.nan_outlier_percent=(s_NaN-s_Mask)/(s-s_Mask)*100;
Cartesian.Mask=CompVel.Mask;
Cartesian.dcor=CompVel.dcor;

end

