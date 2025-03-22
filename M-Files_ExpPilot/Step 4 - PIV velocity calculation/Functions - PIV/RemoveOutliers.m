function [delx] = RemoveOutliers(CompVel)

Threshold=0.5;
delx=CompVel.delta_x;
delz=CompVel.delta_z;
dcor=CompVel.dcor;

delx(dcor < Threshold) = NaN;
delz(dcor < Threshold) = NaN;

s=size(dcor); s=s(1)*s(2);
n=sum(sum((dcor < Threshold)));


B = smoothn(delx, 0.4, 'robust'); % Smooth, remove outliers and replace NaNs
B(isnan(delx)) = NaN; % Put NaNs back

delx(abs(B - delx) > 1) = NaN; % NaN outliers

end

