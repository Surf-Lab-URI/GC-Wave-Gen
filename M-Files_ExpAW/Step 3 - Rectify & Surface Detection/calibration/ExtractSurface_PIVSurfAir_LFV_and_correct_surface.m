function [BadFramePIVSurfLFV,XLFV_Surface,LFV_Surface,XPIV_LFV_Surface,PIV_LFV_Surface,PIV_Surface] = ExtractSurface_PIVSurfAir_LFV_and_correct_surface(PIVSurfA_CamAngle,PIV1_A,idx)

%% PIVSurf Air - LFV CamAngle Correction

%% Step 1: Find surface
% Extract surface
XX = 2001:4000;
YY = 751:7750;%501:8000;
% [imSurf] = Copy_of_FindSurface(PIVSurfA_CamAngle(1001:3000,:), 5, 5);

%%% Pre-extrapolation
if idx > 30000
    PIVSurfA_CamAngle = adapthisteq(PIVSurfA_CamAngle/max(PIVSurfA_CamAngle(:)),'NBins',9);  % increase the contrast and reduce the sensitivity to reduce the noise
end
% When there is a lot of fog, the above pre-processing works well. But not
% always. Check if there are cases when it is convenient to use it.

DX = 75;
S_T = PIVSurfA_CamAngle(XX(1+DX:end),YY);
S_B = PIVSurfA_CamAngle(XX(1:end-DX),YY);
S2 = S_T-S_B;
S = imfilter(S2,fspecial('gaussian',64,64),'replicate');
S(S<0) = NaN;
S = fillmissing(S,'nearest'); %faster than S = smoothn(S,0);

%OG from DE:
% Step = 20;
% Sigma = 5; %20;
% [imSurf] = FindSurface(S, Step, Sigma);

surfSigmas = [100,50,40];
surfSteps = [100,50,];
surfMask = 1;
[imSurf] = CrapperOptimized_FindSurface(S,surfSigmas, surfSteps, surfMask);

% figure(6)
% plot(imSurf.surface)
% imagesc(S)
% colormap gray
% axis equal
% hold on
% plot(imSurf.surface)


PIVSurf_Surface_Raw = imSurf.surface;
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Raw=despike_fab(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = filt_spray(PIVSurf_Surface_Raw);
PIVSurf_Surface_Int = smoothn(PIVSurf_Surface_Int, 'robust');
[SP,~] = spaps(1:length(PIVSurf_Surface_Int), PIVSurf_Surface_Int, 2d3); %Originally 2d4 in DE    %10000);
if length(SP.coefs)>2
    PIVSurf_Surface_A = SP.coefs(2:end-1);
else
    A = polyfit([YY(1) YY(end)], [SP.coefs(1) SP.coefs(end)],1);
    PIVSurf_Surface_A = polyval(A,YY);
end
Usurf = YY;
Vsurf = PIVSurf_Surface_A+XX(1)+DX-1;

%% Step 2 : Resize PIVSurf Air  images to the Size of PIVAir images
U2 = [ 2872 3085 ; 2883 1890 ; 5051 1913 ; 5038 3105]; % The coordinate
% of four points in the PIVSurf Air - LFV surface image (the grid calibration image).
X2 = [ 595 2033 ; 616 296 ; 3781 322; 3768 2062 ]; % The coordinates of  the four
% points in the PIV image. The physical location  of these points are
% the same as PIV surface images. The calibration grid was used to find out
% the exact same locations for these two images.

T2 = maketform('projective',U2,X2);

%%% This is the actual resized image; we retrieve it only for checking,
%%% otherwise we transform only the surface (the resized image is huge!!)
% [Resized_PIVSurf,XPos,YPos] =  imtransform(PIVSurfA_CamAngle,T2,'XYScale',1);
%%% XPos and YPos retrieved from resized transformation!!
XPos(1) = -3.622988898471232e+03;
XPos(2) = 8.903011101528768e+03;
YPos(1) = -2.473308744686341e+03;
YPos(2) = 6.771691255313659e+03;
[Xsurf2,Ysurf2] = tformfwd(T2,Usurf,Vsurf);

% Further rotation to match flat surface
RotAngle = atan((2119.84-2122.86)/(8062.82+2847.82));
M = [cos(RotAngle) sin(RotAngle); -sin(RotAngle) cos(RotAngle)];
Ysurf2 = M*[Xsurf2;Ysurf2];
Xsurf2 = Ysurf2(1,:);
Ysurf2 = Ysurf2(2,:);
Ysurf2 = Ysurf2+17.14;

Xsurf = round(Xsurf2(1):Xsurf2(end));
Xsurf = -2485:7774; % hard-coded to make all the surfaces exactly the same length
Ysurf = interp1(Xsurf2,Ysurf2,Xsurf);

%% Correction of the surface based on the phase
% if idx>10400
%     NBin = 360; % Number of bins in [0 2*pi]
%     DPh = 2*pi/NBin;
%     DMax = 10; % Maximum displacement (corresponding to the crest)
%     XX = 0:DPh:2*pi;
%     %%% Linear displacement (maximum at the crest)
%     % YY = [0:DMax/(NBin/2):DMax-DMax/(NBin/2) DMax DMax-DMax/(NBin/2):-DMax/(NBin/2):0];
%     %%% Gaussian distribution (maximum at the crest)
%     Sigma = 0.25;
%     Mu = pi;
%     C = DMax;
%     YY = C*exp(-0.5*((XX-Mu)/Sigma).^2);
%     % Find Phase
%     Phase2 = angle(hilbert(Ysurf-mean(Ysurf,'omitnan')));
%     % [~,X0(1)] = min(abs(Xsurf-1));
%     % [~,X0(2)] = min(abs(Xsurf-4140));
%     X0 = [1 length(Xsurf)];
%     Phase3 = Phase2(X0(1):X0(2));
%     Phase = wrapTo2Pi(Phase3);
%     Ysurf3 = Ysurf(X0(1):X0(2));
%     Ysurf_CRR = zeros(1,length(Ysurf3));
%     % Apply correction to the surface based on the phase
%     for ii = 1:NBin
%         P = find(Phase>XX(ii) & Phase<XX(ii+1));
%         Ysurf_CRR(P) = Ysurf3(P)-YY(ii);
%     end
%
%     Ysurf = Ysurf_CRR;
%
%     % XPos and YPos are a two-element, real vector that  together  specifiy the
%     % spatial location of the output image B in the 2D output space XY. The two
%     % elements of XPos and YPos give the x-coordinates and y-coordinates of the
%     % first and last columns of B, respectively.
% end
%% Check if bad frame
BadFramePIVSurfLFV = 0;
if imSurf.badFrameBool == 1
    BadFramePIVSurfLFV = 1;
end

%% Results
%%% LFV Surface in LFV coordinates ( == to PIVSurf Air)
XLFV_Surface = Usurf;
LFV_Surface = Vsurf;

%%% LFV Surface in PIV coordinates
XPIV_LFV_Surface = Xsurf;
PIV_LFV_Surface = Ysurf;

%%% PIV Surface from LFV
[~,Ix] = min(abs(Xsurf-1));
[~,Iy] = min(abs((YPos(1):YPos(2))-1));
PIV_Surface = Ysurf(Ix:Ix+size(PIV1_A,2)-1);
%%% Resized PIVSurfA_LFV matching PIV image
% PIV_PIVSurfA = PIVSurfA_LFV_Corrected.img(Iy:Iy+3090-1,Ix:Ix+4140-1);

%% Remind : PIVSurf Air == LFV !!