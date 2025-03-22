function [utSmth, wtSmth, xU] = surface_accel(surfLong,xLong,lambda_cut,params)

%---------------------------------------------------------------------%
% This function accepts a surface wave displacement field and returns %
% the accelerations associated with it...                             %
%                                                                     %
% Call it like:                                                       %
%    utSurf, wtSurf, xAccel = surface_accel(...)                      %
%---------------------------------------------------------------------%

% Unpack parameters
g = params.g;
rho_w = params.WATER_DENSITY;
sigma = params.SURFACE_TENSION;
H = params.WATER_DEPTH;

% Basic properties of the arrays
%x_num = lenght(x_surf)
xL_num = length(xLong);
dx = xLong(2)-xLong(1);

% Butterworth filter to filter out wavelengths shorter than 'lambda_cut':
fs = 1/dx;               % sampling wavenumber [cyc/m]
nyq = fs/2;              % Nyquist wavenumber [cyc/m]
wN = 1/(nyq*lambda_cut); % cutoff wavenumber as fraction of Nyquist
[a,b] = butter(2,wN);    % CHECK a AND b IN CONVERTING TO MATLAB

% Smooth water surface according to the specified Butterworth filter:
surf_smth = filtfilt(a,b,surfLong); % CHECK a AND b IN CONVERTING TO MATLAB
% % % surf_smth = surfLong;

% Use FFT to approximate water surface and its velocity: ------------
% Compute wave numbers in fft:
fft_freq = 0:fs/xL_num:nyq;
% % % fft_freq = 0:fs/xL_num:fs-fs/xL_num;
k = 2*pi*fft_freq';        % wavenumber in [rad/m]

% Use dispersion relationship to compute angular frequency in [rad/s]:
omega = sqrt(k.*(g + sigma*k.^2/rho_w).*tanh(k*H));

% Adjustment for capillary waves
k_cut = 2*pi/0.004;
omega(k > k_cut) = k(k > k_cut).*sqrt(sigma/rho_w*k(k > k_cut));

% Do fast fourier transform:
Fft = fft(surf_smth);
Fft(2:length(Fft)/2) = 2*Fft(2:length(Fft)/2);
Fft(length(Fft)/2+2:end) = [];

% Extract sine and cosine coefficients from the fft:
cosCoef = real(Fft)';
sinCoef = -imag(Fft)';

%compute sine and cosine coefficients of orbital velocity component u
u_cosCoef = (cosCoef.*omega).*tanh(k*H);      % CHECK consistency with python
u_sinCoef = (sinCoef.*omega).*tanh(k*H);      % CHECK consistency with python

%compute sine and cosine coefficients of orbital velocity component w
w_cosCoef = -sinCoef.*omega;      % CHECK consistency with python
w_sinCoef = cosCoef.*omega;      % CHECK consistency with python

% inverse fft in order to get orbital velocity components
u = real(ifft(u_cosCoef - 1i*u_sinCoef, xL_num));
w = real(ifft(w_cosCoef - 1i*w_sinCoef, xL_num));

% remove noisy region on side
% re-apply low-pass filter
chop = 200;
xU = xLong(chop+1:end-chop);
u = filtfilt(a, b, u(chop+1:end-chop));
w = filtfilt(a, b, w(chop+1:end-chop));
% % % u = u(chop+1:end-chop);
% % % w = w(chop+1:end-chop);
xU_num = length(xU);

% Use FFT on orbital velocities to get surface acceleration -------------------
%Compute wave numbers in fft for shortened domain
fft_freq = 0:fs/xU_num:nyq;
% % % fft_freq = 0:fs/xU_num:fs-fs/xU_num;
k = 2*pi*fft_freq';        % wavenumber in [rad/m]

% Use dispersion relationship to compute angular frequency [rad/s]:
omega = sqrt(k.*(g + sigma*k.^2/rho_w).*tanh(k*H)) ;

%Adjustment for capilary waves
k_cut = 2*pi/0.004;
omega(k > k_cut) = k(k > k_cut).*sqrt(sigma/rho_w*k(k > k_cut));

% compute fft of u and w orbital components
fftU = fft(u);
fftU(2:length(fftU)/2) = 2*fftU(2:length(fftU)/2);
fftU(length(fftU)/2+2:end) = [];
fftW = fft(w);
fftW(2:length(fftW)/2) = 2*fftW(2:length(fftW)/2);
fftW(length(fftW)/2+2:end) = [];

% extract sine and cosine components
uSmth_cosCoef = real(fftU);
uSmth_sinCoef = -imag(fftU);

wSmth_cosCoef = real(fftW);
wSmth_sinCoef = -imag(fftW);

% compute sine and cosine componets of the accelerations
utSmth_cosCoef = -uSmth_sinCoef.*omega;
utSmth_sinCoef = uSmth_cosCoef.*omega;

wtSmth_cosCoef = -wSmth_sinCoef.*omega;
wtSmth_sinCoef = wSmth_cosCoef.*omega;


% Invert FFT to get surface accelrations
ut = real(ifft(utSmth_cosCoef - 1i*utSmth_sinCoef, xU_num));
wt = real(ifft(wtSmth_cosCoef - 1i*wtSmth_sinCoef, xU_num));

% smooth surface accelerations
utSmth = filtfilt(a, b, real(ut));
wtSmth = filtfilt(a, b, real(wt));
% pause

%------------------------------ End function: ------------------------------%
%---------------------------------------------------------------------------%