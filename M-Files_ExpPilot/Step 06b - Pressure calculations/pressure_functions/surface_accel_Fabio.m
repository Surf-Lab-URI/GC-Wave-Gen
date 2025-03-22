function [utSmth, wtSmth, xU] = surface_accel_Fabio(surfLong,xLong,lambda_cut,params)

%---------------------------------------------------------------------%
% Compute acceleration vector from free surface elevation             %
%   utSmth = smoothed acceleration in x                               %
%   wtSmth = smoothed acceleration in z                               %
%                                                                     %
% Call it like:                                                       %
%   [utSurf, wtSurf, xU] = surface_accel_Fabio(...)                   %
%---------------------------------------------------------------------%

% Unpack parameters
g = params.g;
rho_w = params.WATER_DENSITY;
sigma = params.SURFACE_TENSION;
H = params.WATER_DEPTH;

% Basic properties of the arrays
xL_num = length(xLong);
dx = xLong(2)-xLong(1);

% Butterworth filter to filter out wavelengths shorter than 'lambda_cut':
fs = 1/dx;               % sampling wavenumber [cyc/m]
nyq = fs/2;              % Nyquist wavenumber [cyc/m]
wN = 1/(nyq*lambda_cut); % cutoff wavenumber as fraction of Nyquist
[a,b] = butter(2,wN);    

% Smooth water surface according to the specified Butterworth filter:
surf_smth = filtfilt(a,b,surfLong); 
% surf_smth = filtfilt(a,b,surf_smth); % double smooth if surface is very discontinuous 

% Use FFT to approximate water surface and its velocity: ------------
% Compute wave numbers in fft:
fft_freq = 0:fs/xL_num:fs-fs/xL_num;
k = 2*pi*fft_freq;        % wavenumber in [rad/m]

% Use dispersion relationship to compute angular frequency in [rad/s]:
om = sqrt(k.*(g + sigma*k.^2/rho_w).*tanh(k*H));

% Do fast fourier transform:
f = fft(surf_smth);
Sym_Vec = NaN(size(f));
Sym_Vec(1) = 1; Sym_Vec(length(f)/2+1) = 1;
Sym_Vec(2:length(f)/2) = 2;
Sym_Vec(length(f)/2+2:end) = 0;
F = (f.*Sym_Vec);

% Invert FFT to get surface accelerations
surf_acc = filtfilt(a,b,ifft(-F.*om.^2,xL_num,'nonsymmetric'));

% Crop sides
chop = 200;
xU = xLong(chop+1:end-chop)';
ut = -imag(surf_acc(chop+1:end-chop))'; % x component
wt = real(surf_acc(chop+1:end-chop))'; % z component

% smooth surface accelerations
utSmth = filtfilt(a, b, ut);
wtSmth = filtfilt(a, b, wt);

% pause


%------------------------------ End function: ------------------------------%
%---------------------------------------------------------------------------%