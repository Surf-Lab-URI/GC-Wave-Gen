function [Surface] = Surface_from_PIVSURF_to_PIV_h8cm(x)

%% Surface correction to obtain PIV surface with PIVSURF
load Surface_match DeltaS Filt Filt2 MM Alpha M

X = [(1:length(x));x];
Xr = M*X;
xr = interp1(Xr(1,:),Xr(2,:),1:length(x),'pchip','extrap');

Surface = Filt.*xr-MM;