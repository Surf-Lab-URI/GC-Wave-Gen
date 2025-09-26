function wave_numb=wav_num(t,x1)

% N=2^nextpow2(length(t));
% yt1=2/N*abs(fft(x1,N));    % ??N??2?????
% f = Fs/N*(0:1:N-1);
% % ft1=(1:length(t))/length(t)/(t(2)-t(1));
% ft1=f(2:N/2);yt1=yt1(2:N/2)*2;
% semilogx(ft1,yt1)
% grid on
% %xlabel('Frequency/Wave number (1/unit of x)','fontsize',10)
% xlabel('wave number','fontsize',10)
% ylabel('Power','fontsize',10)


% yt1=2/length(x1)*abs(fft(x1));
% 
% ft1=(t(end)-t(1))/length(t);
% ft2=t(end)/2;
% ft=[ft1:(ft2-ft1)/length(t)*2:ft2];
% ft=ft(2:end);
% yt1=yt1(2:(end+1)/2)*2;
% plot(2*pi./ft,yt1)

% yt1=2/length(x1)*abs(fft(x1));
% ft1=(1:length(t))/length(t)/(t(2)-t(1));
% ft1=ft1(2:end/2);yt1=yt1(2:end/2)*2;
% % semilogx(1./ft1,yt1);
% [y_max index]=max(yt1);
% wave_numb=2*pi*ft1(index);

sp=spectrum(x1,length(x1),0,hanning(length(x1)));
sp=sp(:,1);
 DX=5.65571e-05;
k=[0:length(sp)-1]/(length(sp)-1)*(2*pi/2/DX);
[y_max index]=max(sp);
wave_numb=k(index);
end

% 
% Fs = 1000;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% L = 1000;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
% X = S + 2*randn(size(t));
% plot(1000*t(1:50),X(1:50))
% title('Signal Corrupted with Zero-Mean Random Noise')
% xlabel('t (milliseconds)')
% ylabel('X(t)')
% 
% Y = fft(X);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% plot(f,P1)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')

