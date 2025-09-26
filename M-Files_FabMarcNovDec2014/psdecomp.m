function [Pu,Pv,Su,Sv] = psdecomp(uu,vv);
% Separate (uu,vv) into irrotational (potential flow) and
% nondivergent (Solenoidal flow) parts. This implementation
% assumes dx = dy, though the arrays can be rectangular.
% original implementation by J. A. Smith, May 23, 2005;
% function updated by J. A. Smith, June 28, 2007.
[mn, nn] = size(uu); % vv is the same size.
nft = 2^nextpow2(max(mn,nn)*1.2);% min 20% zero pad.
nf2 = nft/2;
kv = (-nf2:(nf2-1))/nft;
thetamn= fftshift(atan2(kv(ones(nft,1),:)',kv(ones(nft,1),:)));
vk1 = fft2(uu,nft,nft);
vk2 = fft2(vv,nft,nft);
vkdiv = vk1.*cos(thetamn)+vk2.*sin(thetamn);
vav = real(ifft2(vkdiv.*cos(thetamn)));
Pu = vav(1:mn,1:nn);
vav = real(ifft2(vkdiv.*sin(thetamn)));
Pv = vav(1:mn,1:nn);
vk1 = fft2(uu,nft,nft);
vk2 = fft2(vv,nft,nft);
vkdiv = vk2.*cos(thetamn)-vk1.*sin(thetamn);
vav = real(ifft2(-vkdiv.*sin(thetamn)));
Su = vav(1:mn,1:nn);
vav = real(ifft2(vkdiv.*cos(thetamn)));
Sv = vav(1:mn,1:nn);