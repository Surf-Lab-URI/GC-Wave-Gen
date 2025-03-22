function [Cartesian] = RemoveOutliers_fab(CompVel)

Threshold=0.4;
delx=CompVel.delta_x;
delz=CompVel.delta_z;
dcor=CompVel.dcor;

delx(dcor < Threshold) = NaN;
delz(dcor < Threshold) = NaN;

s=size(dcor); s=s(1)*s(2);
nc=sum(sum((dcor < Threshold)));
ThrX = (delx-CompVel.INTdelx)>4 | (delx-CompVel.INTdelx)<-6;
ThrZ = (delz-CompVel.INTdelz)>2 | (delz-CompVel.INTdelz)<-2;

nx=sum(sum(abs(delx-CompVel.INTdelx)>1));
nz=sum(sum(abs(delz-CompVel.INTdelz)>1));

delx(ThrX) = NaN; % NaN outliers
delz(ThrZ) = NaN; % NaN outliers



%% Outlier removed
Cartesian.delx=delx.*CompVel.Mask;
Cartesian.delz=delz.*CompVel.Mask;

%% Smoothed
Cartesian.u=smoothn(Cartesian.delx,0.4,'robust');
Cartesian.w=smoothn(Cartesian.delz,0.4,'robust');
Cartesian.corr_outlier_percent=nc/s*100;
Cartesian.x_outlier_percent=nx/s*100;
Cartesian.z_outlier_percent=nz/s*100;
Cartesian.Mask=CompVel.Mask;
Cartesian.dcor=CompVel.dcor;

end

