function [xWG2,heightWG2,xWG4,heightWG4] = extract_WG_signal(imagenameWG2,imagenameWG4)

%%


%%

% Parameters for image loading
nl = 472; % number of columns used
Sigma = 5; %sigma for gaussian filter
Size = 5; %kernel size for gaussina filter
fileoffset = 0;  %file header length in bytes using Dynamic Studio

%% WG12
%%% For shoaling waves, only WG2 used (WG1 removed)

%%% WG12 image from :
%\\spray1\d\Shoaling_waves\data\Calibration\05082022_calibration\Movie1_Scene1\RAW\WG12\
%%% PIVSURF image from :
%\\spray1\d\Shoaling_waves\data\Calibration\05082022_calibration\Movie2_Scene2\RAW\PIVSURF\Movie2_Scene2_PIVSURF_06.raw

%%% Load image

fid = fopen( imagenameWG2);
fseek( fid , fileoffset , 'bof' );
Img = fread(fid,nl*2048,'uint16');
Img = rot90(reshape(Img,2048,nl),2);
fclose(fid);

% Make Img2 homogeneous in the whole image 
Img2 = Img;
Img2(1:1024,:) = Img(1:1024,:)-mean(mean(Img(1:1024,:)));
Img2(1025:end,:) = Img(1025:end,:)-mean(mean(Img(1025:end,:)));

% Filter out noise
img = imgaussian(Img2,Sigma,Size);
% figure;imagesc(img);colormap gray
% axis equal


%%% WG12 Transformation

% Transform WG12 image in PIV coordinates
U0 = fliplr([ 1239 107 ; 1164 107 ; 1164 32 ; 1239 32 ]); % Points in WG12 coordinates
X0 = fliplr([ 5753 2041 ; 5327 2040 ; 5325 1613 ; 5751 1615]); % Points in PIV coordinates

%Dpix_WG12 = mean(nanmean(dist(X0')*59.476d-6./dist(U0'))); % WG12 resolution (in m/pix)

TForm = maketform('projective',U0,X0);

[Resized_WG12,XPos,YPos] =  imtransform(img,TForm,'XYScale',1);
X_WG12 = (1:size(Resized_WG12,2))+XPos(1)-2468;
Y_WG12 = (1:size(Resized_WG12,1))+YPos(1);

%imagesc(X_WG12,Y_WG12,img)

%%% Obtain single-point surface for WG2

X2 = 0;
Int = 701:1200; % 801:950;
smth_vert= smoothn(Resized_WG12(:,X2+1:end),10000);
gv = gradient(mean(smth_vert(:,Int),2,'omitnan'));
% [~,eta_locs2] = findpeaks(gv,'npeaks',1,'SortStr','descend');

% Find maximum around the smoothed version of smth_vert to avoid noise
% We also know that eta_locs2_ must be ~8000; so, if it's too small or too big, it's wrong
gv2 = gradient(smooth(mean(smth_vert(:,Int),2,'omitnan'),150));
eta_locs2_ = 2;
while eta_locs2_<5000 
    gv2(1:eta_locs2_) = mean(gv2,'omitnan');
    [~,eta_locs2_] = findpeaks(gv2,'npeaks',1,'SortStr','descend');
    if eta_locs2_>10000
        gv2(eta_locs2_:end) = mean(gv2,'omitnan');
        eta_locs2_ = 2;
    end
end

% Now retrieve the exact maximum of gv around gv2
% [~,eta_locs2] = max(gv(eta_locs2_-100+1:eta_locs2_+100));
[~,eta_locs2] = findpeaks(gv(eta_locs2_-200+1:eta_locs2_+200),'minpeakheight', max(gv(eta_locs2_-200+1:eta_locs2_+200))/2, 'npeaks',1);
if isempty(eta_locs2)
    eta_locs2 = eta_locs2_;
else
    eta_locs2 = eta_locs2+(eta_locs2_-200);
end

[~,x_locs2] = max(mean(Resized_WG12(eta_locs2:eta_locs2+500,:)));
% [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/2, 'npeaks',1);

heightWG2 = Y_WG12(eta_locs2);
xWG2 = X_WG12(x_locs2);

%% WG34
%%% For shoaling waves, only WG4 useful (WG3 intensity very low, difficult to separate from the background noise)

%%% Load image
n2 = 159;

fid = fopen( imagenameWG4);
fseek( fid , fileoffset , 'bof' );
Img = fread(fid,nl*2048,'uint16');
Img = reshape(Img,2048,nl);
fclose(fid);

% Make Img2 homogeneous in the whole image 
Img2 = Img;
Img2(1:1024,:) = Img(1:1024,:)-mean(mean(Img(1:1024,:)));
Img2(1025:end,:) = Img(1025:end,:)-mean(mean(Img(1025:end,:)));

% Filter out noise
img = imgaussian(Img2,Sigma,Size);
% figure;imagesc(img);colormap gray
% axis equal

%%% WG34 Transformation

% Transform WG34 image in PIV coordinates
U0 = [ 80 793 ; 156 794 ; 156 717 ; 80 717 ]; % Points in WG34 coordinates
X0 = [ 8444 4927 ; 8869 4927 ; 8874 4499 ; 8447 4498]; % Points in PIV coordinates

%Dpix_WG34 = mean(nanmean(dist(X0')*59.476d-6./dist(U0'))); % WG34 resolution (in m/pix)

TForm = maketform('projective',U0,X0);

[Resized_WG34,XPos,YPos] =  imtransform(img,TForm,'XYScale',1);
DX = -2468;
X_WG34 = (1:size(Resized_WG34,2))+XPos(1)+DX;
Y_WG34 = (1:size(Resized_WG34,1))+YPos(1);

%%% Obtain single-point surface for WG4

X4 = 850;
Int = 100:300;
smth_vert= smoothn(Resized_WG34(:,X4+1:end),10000);
gv = gradient(mean(smth_vert(:,Int),2,'omitnan'));
% [~,eta_locs4] = findpeaks(gv,'npeaks',1,'SortStr','descend');

% Find maximum around the smoothed version of smth_vert
gv2 = gradient(smooth(mean(smth_vert(:,Int),2,'omitnan'),100));
eta_locs4_ = 2;
while eta_locs4_<3000 
    gv2(1:eta_locs4_) = mean(gv2,'omitnan');
    [~,eta_locs4_] = findpeaks(gv2,'npeaks',1,'SortStr','descend');
    if eta_locs4_>10000
        gv2(eta_locs4_:end) = mean(gv2,'omitnan');
        eta_locs4_ = 2;
    end
end

% Now retrieve the exact maximum of gv around gv2
% [~,eta_locs4] = max(gv(eta_locs4_-100+1:eta_locs4_+100));
[~,eta_locs4] = findpeaks(gv(eta_locs4_-200+1:eta_locs4_+200),'minpeakheight', max(gv(eta_locs4_-200+1:eta_locs4_+200))/2, 'npeaks',1);
if isempty(eta_locs4)
    eta_locs4 = eta_locs4_;
else
    eta_locs4 = eta_locs4+(eta_locs4_-200);
end

[~,x_locs4] = max(mean(Resized_WG34(eta_locs4:eta_locs4+500,X4+1:end)));
% [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/2, 'npeaks',1);

heightWG4 = Y_WG34(eta_locs4);
xWG4 = X_WG34(x_locs4+X4);

