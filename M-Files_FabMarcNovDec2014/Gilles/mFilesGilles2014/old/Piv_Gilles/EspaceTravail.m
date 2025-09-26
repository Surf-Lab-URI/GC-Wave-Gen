
clear all
close all

%% Chargement de l'image

% Image eau calme
% load('F:\data\ExpLC1_dt25ms_1\rawImages\Lfv\ExpLC1_dt25ms_1_Lfv_0001.mat')
% Surface environ au pixel 1154
% Image au temps longs
load('F:\data\ExpLC1_dt25ms_1\rawImages\Lfv\ExpLC1_dt25ms_1_Lfv_2151.mat')

%% Detection de surface
% Methode Gilles
tic,
img=imgLfv;
surfEdge=edge(img,'sobel');
surfEdgeSmt=5*smoothn(surfEdge);
for i=1:size(surfEdgeSmt,2)
    j=1;
%     if surfEdgeSmt(j
% if isempty(find(surfEdgeSmt>=0.9))
%     surface(i) = nan;
%     return
% end
    while surfEdgeSmt(j,i)<0.9
        j=j+1;
        if j>2048
            surface(i) = nan;
        end
    end
    surface(i)=j;
end
surface=smoothn(surface,100);
figure,imagesc(img), colormap(bone)
hold on, plot(surface,'g')
smth_surface = surface
toc

% Methode Marc

disp('Recherche de la surface (Marc)');
tic,
img = imgLfv;
smth_vert = nan(size(img));
grad_vert = nan(size(img));
mask = nan(size(img));
surface = nan(1,size(img,2));
%
for i=1:size(img,2)
%     keyboard
    imgi = img(:,i);
    %locate outliers and nan them
    imgistd = std(imgi);
    imgimean = mean(imgi);
    imgi(abs(imgi-imgimean)>3*imgistd)=nan;
    %smooth each column
    smth_vert(:,i) = smoothn(imgi,10000);
    %compute gradient on each column
    grad_vert(:,i) = gradient(smth_vert(:,i));
    gv = grad_vert(:,i);
    [~,locs] = findpeaks(gv,'minpeakheight', max(gv)/2, 'npeaks',1);
    surface(i) = locs;
    clear imgi gv locs pks
    clear grad_vert %%%%%%%%%%%%%%%%%%%%%%%%
end
%
%locate outliers and nan them
surfstd = std(surface);
surfmean = mean(surface);
surface(abs(surface-surfmean)>3*surfstd)=nan;
smth_surface = smoothn(surface);
clear surfstd surfmean surface smth_vert %%%%%%%%%%%%%%%
toc
%figure,imagesc(img), colormap(bone)
% hold on, plot(surface,'r')

%% Comparaison surface image corrige et non corrige
surf11 = surf1*pivsurf_res/piv_res;
pivsurf_res = 102e-6; %um/pixel
piv_res = 40e-6; %um/pixel
x = 0:pivsurf_res:2047*pivsurf_res;  % initial data sites
xx = 0:piv_res:3784*piv_res; % target data sites
% ssrvs = spline(x,ssrv,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
surf1_interp = spline(x,surf11,xx);  % smth_surface_rescaled_vert_splined (spline interpolation)
surf1_fin = surf1_interp - 2226;
figure, plot(surf1_fin(847:end)), hold on, plot(surf2, 'r')


%% Ajout d'un mask
for i=1:size(img,2)
    mask(round(smth_surface(i):2048),i) = 1;
end

%% Premiere estimation de la vitesse (calcul des vitesses orbitales grace a la decomposition em mode de la surface de la Lfv)

disp('Calcul des orbitales')
tic,
s=smth_surface;
dx=1;


x_s=[0:8:2048]*dx;
z_s=[0:8:1024]*dx;



u=zeros(length(z_s)-1,length(x_s)-1);w=zeros(length(z_s)-1,length(x_s)-1);
f=fft(s); %Fourrier Modal decomposition
fa=2*abs(f)/length(f); %amplitude of mode
fp=angle(f); %phase of mode
k=[0:(length(s)-1)]/(length(s)-1)*2*pi/dx; %wavenumber of mode
z_sc=z_s-0; %rectification niveau surface
for j=1:length(z_s)-1
g=0;
h=0;
for i=1:floor(length(s)/2)%/10
g=g+fa(i)*sqrt(9.81*k(i))*exp(-z_sc(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i))-mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_sc(j)*k(i));
h=h+fa(i)*sqrt(9.81*k(i))*exp(-z_sc(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i))-mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_sc(j)*k(i));
end
u(j,:)=g;
w(j,:)=h;
end

for i=1:length(s) %On place pour chaque colonne le z=0 au niveau de la surface
    u(:,i)= 
end



% % Version initiale
% x_s=[0:8:2048];
% z_s=[0:8:2048];
% 
% u=zeros(length(z_s),length(x_s)-1);w=zeros(length(z_s)-1,length(x_s)-1);
% f=fft(s); %Fourrier Modal decomposition
% fa=2*abs(f)/length(f); %amplitude of mode
% fp=angle(f); %phase of mode
% k=[0:(length(s)-1)]/(length(s)-1)*2*pi/dx; %wavenumber of mode
% for j=1:length(z_s)
% g=0;
% h=0;
% for i=1:floor(length(s)/2)%/10
% g=g-fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i))+mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
% h=h+fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i))-mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
% end
% u(j,:)=g;
% w(j,:)=h;
% end
% 






figure, imagesc(mask)
hold on,  quiver(10:10:2048,10:10:2028,u(10:10:2028,10:10:2048),w(10:10:2028,10:10:2048))


dim1=256
dim2=256
figure, quiver(2:dim2,2:dim1,u(2:dim1,2:dim2),w(2:dim1,2:dim2))


% % Visualisation vorticite

Vort=curl(u,w);
figure, imagesc(Vort), caxis([0,0.01])

%%







%% Surface RECONSTRUCTION at t
x_s=[1:2048];
dx=1;
t=1/7.2;
f=fft(s);
fa=2*abs(f)/length(f);
fp=angle(f);
k=[0:(length(s)-1)]/(length(s)-1)*2*pi/dx;
omega=sqrt(9.81*k);
g=0;
for i=1:floor(length(s)/2)
g=g+(fa(i))*cos(k(i)*x_s(1:end-1)+fp(i)-omega(i)*t)-mean(s)/floor(length(s)/2);
end
g=interp1(x_s(1:end-1),g,x_s,'linear','extrap');%X_s(2:end)

figure, plot(s)
hold on, plot(g,'r')










%
%   g=isnan(mean(stack_of_profiles'));
%  b=ones(1,ceil(1.7d-2*1/dx))/ceil(1.7d-2*1/dx);
% CC=[];
% PP=[];
%
%  for i=1:J-1
%      if (g(i)+g(i+1))<1 %two in a row without bad profiles
%
%          %PHASE SPEED
%          s1=delx(i,:); lowpass_s1=filtfilt(b,1,s1); % take a profile then the next one
%          s2=delx(i+1,:);lowpass_s2=filtfilt(b,1,s2);
%          [c_cov,lags_cov] = xcov(lowpass_s1,lowpass_s2,'coeff');
%          c_max = max(c_cov(1:length(s1)));
%          lag_off = -lags_cov(find(c_cov==c_max,1,'first'));
%          c=dx*lag_off/dt;
%          C(i)=c;
%          quality_of_phase(i)=c_max;
%          [~,locs] = findpeaks(lowpass_s1);
%          a_crest(i)=mean(s1(locs));
%          a_crest_error(i)=std(s1(locs));
%
%
%
% %WAVELENGTH & AMPLITUDE
% if length(locs)>1
%     lambda_crest(i)=mean(diff(locs))*dx;
%     lambda_crest_error(i)=std(diff(locs))*dx;
% else
%     lambda_crest(i)=NaN;
%     lambda_crest_error(i)=NaN;
% end
% [~,locs] = findpeaks(-lowpass_s1);
%     a_trough(i)=mean(s1(locs));
%     a_trough_error(i)=std(s1(locs));
% if length(locs)>1
%     lambda_trough(i)=mean(diff(locs))*dx;
%     lambda_trough_error(i)=std(diff(locs))*dx;
% else
%     lambda_trough(i)=NaN;
%     lambda_trough_error(i)=NaN;
% end
% z=find(abs(diff(sign(lowpass_s1)))>1);
% if length(z)>1
%     lambda_zero(i)=mean(diff(z))*dx*2;
%     lambda_zero_error(i)=std(diff(z))*dx*2;
% else
%     lambda_zero(i)=NaN;
%     lambda_zero_error(i)=NaN;
% end
%
%      else
%         C(i)=NaN;
%         quality_of_phase(i)=NaN;
%         lambda_crest(i)=NaN;
%         lambda_crest_error(i)=NaN;
%         a_crest(i)=NaN;
%         a_crest_error(i)=NaN;
%         lambda_trough(i)=NaN;
%         lambda_trough_error(i)=NaN;
%         a_trough(i)=NaN;
%         a_trough_error(i)=NaN;
%         lambda_zero(i)=NaN;
%         lambda_zero_error(i)=NaN;
%
%      end
%  end
%
%
%
%
%
%
