function [Surface] = Reconstruct_Surface_h8cm(x,y)

%% Surface correction to obtain PIV surface with PIVSURF
Dm = 26.6294;
% Rotate x
cc = polyfit(1:size(x,2),x-y-Dm,1);
xr = x;
while abs(cc(2))>0.1
    X = [(1:length(xr));xr];
    Alpha = rad2deg(-atan(cc(1)));
    M = [cosd(Alpha),-sind(Alpha);sind(Alpha),cosd(Alpha)];
    Xr = M*X;
    xr = interp1(Xr(1,:),Xr(2,:),1:size(x,2),'pchip','extrap');
    cc = polyfit(1:size(x,2),xr-y-mean(xr-y,'omitnan'),1);
end
%figure; plot(xr-y-nanmean(xr-y),'b')
Dmr = 31.9781;

% Define filter1
Filt1(1:340) = 1.3;
Filt1(341:514) = (1.3:-(1.3-1.1323)/(514-341):1.1323);
Filt1(515) = 1.0425;
Filt1(516:800) = 1.0425+(1.08-1.0425)/(800-515):(1.08-1.0425)/(800-515):1.08;
Filt1(801:996) = 1.06-(1.06-1.0)/(996-800):-(1.06-1.0)/(996-800):1.0;
Filt1(781:820) = Filt1(780)+(Filt1(821)-Filt1(780))/(821-780):(Filt1(821)-Filt1(780))/(821-780):Filt1(820);

% Define filter2
XX = xr-mean(xr-y);XX= XX(1201:2100);
YY = y; YY = YY(1201:2100);
Filt2(1:900) = XX./YY;

% Reconstruct surface
xIn = 2411;
xFin = 3407;
Surface([1:xIn-1,xFin:length(xr)]) = xr([1:xIn-1,xFin:length(xr)])-Dmr;
Surface(2411:3406) = interp1(Filt1.*(1:996)+(xIn-1),xr(2411:3406)-Dmr,(1:996)+(xIn-1),'pchip','extrap');
Surface(1201:2100) = Filt2.*Surface(1201:2100);
%figure;plot(Surface,'r');hold on;plot(y,'b')