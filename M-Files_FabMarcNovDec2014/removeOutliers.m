function delx = removeOutliers( u, dcor )
delx = u;
delx(dcor < 0.5) = NaN;
%
B = smoothn(delx,0.4, 'robust'); % smoothe, remove outliers and replace nans
B(isnan(delx)) = nan; % put nans back
%
delx(abs(B-delx)>1) = nan; % nan outliers
end