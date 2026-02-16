

d1=imresize(imgPivsurf1,176.77/103.48); %Resizing PIVsurf image to match PIV image (1/17677 m per pix for PIV and 1/10348 m per pix for PIVsurf)
%d1 is PIVsurf image with the same resolution as the PIV image)
Surface_PIV=Surface_PIV1(724:724+2047); %in this resized surface image, x=724 is the left edge of PIV image
%Surface_PIV1 is the surface detected on teh PIVsurf image and 724:724+2047 is the portion that corresponds to the PIV image
Surface_PIV=Surface_PIV1-1838+370; %still water is @ Y=1838 in PIVsurf and Y=370 in PIV - This locates the surface in the PIV image so I can make a mask


dx =1/17677; %m per pix
dt=10d-3; %second per pair
gravity=9.81; %m/s^2
gravity_piv=gravity/dx*dt^2; %gravity in pix per pair


f = fft(Surface_PIV1); 
fa = 2*abs(f)/length(f); %amplitude of mode
fp = angle(f);
k = [0:(length(Surface_PIV1)-1)]/(length(Surface_PIV1)-1)*2*pi;

z_s = 0:1:2048-1; %depth at which we'll estimate orbital vel
x_s = 724:1:724+2047; %x at which we'll estimate orbital vel
for j = 1:length(z_s)
    g = 0;
    h = 0;
    su=0;
    for i = 1:1:floor(length(k)/2/10)% Choose how many mode to keep; here, I keep 10 times less than the max number of modes (i.e. low pass...)
        su=su+fa(i)*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); 
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only low modes for orbital vel from lfv_surf
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
    end
    SU(j,:)=su-z_s(j); 
    u(j,:) = -g;   %minus sign here because the surface is technically "upside down" with low pixvalues at crest - that's because the origin (0,0) of an image is in the upper left
    w(j,:) = h;
end


%Calculating orbital velocity on the grid of the final PIV results
clear u; clear w
GS=compVel.GS;
x_s = 724+GS:GS:724+2047;
z_s = 0:GS:2048-GS*2;
for j = 1:length(z_s)
    g = 0;
    h = 0;
       for i = 1:floor(length(k)/2/GS/10)%/10
        g = g+fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end)+fp(i)); % only low modes for orbital vel from lfv_surf
        h = h-fa(i)*sqrt(gravity_piv*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end)+fp(i));
       end
    u(j,:) = -g;
    w(j,:) = h;
end
