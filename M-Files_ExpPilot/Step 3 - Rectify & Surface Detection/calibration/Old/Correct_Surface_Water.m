function [XPIVWater_Surface,PIVWater_Surface,XPIV_PIVWater_Surface_CRR,PIV_PIVWater_Surface_CRR] = Correct_Surface_Water(PIVWater_CRR,XPIV_PIVSurfW_Surface,PIV_PIVSurfW_Surface,PIV_W,idx)

%% Correction of the surface based on the phase
if idx>400
    NBin = 360; % Number of bins in [0 2*pi]
    DPh = 2*pi/NBin;
    DMax = 50; % Maximum displacement (corresponding to the crest)
    XX = 0:DPh:2*pi;
    %%% Linear displacement (maximum at the crest)
    % YY = [0:DMax/(NBin/2):DMax-DMax/(NBin/2) DMax DMax-DMax/(NBin/2):-DMax/(NBin/2):0];
    %%% Gaussian distribution (maximum at the crest)
    Sigma = 0.75;
    Mu = pi;
    C = DMax;
    YY = C*exp(-0.5*((XX-Mu)/Sigma).^2);
    % Find Phase
    Phase2 = angle(hilbert(PIV_PIVSurfW_Surface-mean(PIV_PIVSurfW_Surface,'omitnan')));
    [~,X0(1)] = min(abs(XPIV_PIVSurfW_Surface-PIVWater_CRR.Xpos(1)));
    [~,X0(2)] = min(abs(XPIV_PIVSurfW_Surface-PIVWater_CRR.Xpos(2)));
    Phase3 = Phase2(X0(1):X0(2));
    Phase = wrapTo2Pi(Phase3);
    Surface_Water = PIV_PIVSurfW_Surface(X0(1):X0(2));
    Surface_Water_CRR = zeros(1,length(Surface_Water));
    % Apply correction to the surface based on the phase
    for ii = 1:NBin
        P = find(Phase>XX(ii) & Phase<XX(ii+1));
        Surface_Water_CRR(P) = Surface_Water(P)+YY(ii);
    end
else 
    [~,X0(1)] = min(abs(XPIV_PIVSurfW_Surface-PIVWater_CRR.Xpos(1)));
    [~,X0(2)] = min(abs(XPIV_PIVSurfW_Surface-PIVWater_CRR.Xpos(2)));
    Surface_Water_CRR = PIV_PIVSurfW_Surface(X0(1):X0(2));
end
% % Check surface position
% figure;imagesc(PIVWater_CRR.Xpos(1):PIVWater_CRR.Xpos(2),PIVWater_CRR.Ypos(1):PIVWater_CRR.Ypos(2),PIVWater_CRR.img);colormap gray
% hold on;plot(XPIV_PIVSurfW_Surface,PIV_PIVSurfW_Surface,'g')
% X = XPIV_PIVSurfW_Surface;
% plot(X(X0(1)):X(X0(2)),Surface_Water_CRR,'r')

%% Further rotation to match PIVWater Surface
RotAngle = -atan(6/4500); %20/4500
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf2 = M*[1:length(Surface_Water_CRR);Surface_Water_CRR];
Xsurf2 = Ysurf2(1,:);
Ysurf2 = Ysurf2(2,:)+10;
XPIV_PIVWater_Surface_CRR = round(Xsurf2(1)):round(Xsurf2(end));
PIV_PIVWater_Surface_CRR = interp1(Xsurf2,Ysurf2,XPIV_PIVWater_Surface_CRR);

%% Surface Water in PIVWater coordinates
RotAngle = -atan(27/4500);
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf = M*[1:length(Surface_Water_CRR);Surface_Water_CRR];
Xsurf = Ysurf(1,:);
Ysurf = Ysurf(2,:)+18;
X = XPIV_PIVSurfW_Surface;
[U2,V2] = tforminv(PIVWater_CRR.Tform,X(X0(1)):X(X0(2)),medmob2(Ysurf',20)');
[~,Ix1] = min(abs(U2-1));
[~,IxEnd] = min(abs(U2-size(PIV_W,2)));
% XPIVWater_Surface = round(U2(1)):round(U2(end));
% PIVWater_Surface = interp1(U2,V2,XPIVWater_Surface,'linear','extrap');
XPIVWater_Surface = round(U2(Ix1)):round(U2(IxEnd));
PIVWater_Surface = interp1(U2,V2,XPIVWater_Surface,'linear','extrap');