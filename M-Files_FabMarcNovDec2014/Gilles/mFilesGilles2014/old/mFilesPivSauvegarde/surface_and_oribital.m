%Surface RECONSTRUCTION

f=fft(s);
fa=2*abs(f)/length(f);
fp=angle(f);
k=[0:(length(s)-1)]/(length(s)-1)*2*pi/dx;
g=0;
for i=1:floor(length(s)/2)
g=g+(fa(i))*cos(k(i)*x_s(1:end-1)+fp(i))-mean(s)/floor(length(s)/2);
end
g=interp1(x_s(1:end-1),g,x_s,'linear','extrap');%X_s(2:end)
g=g-mean(s);

%Surface RECONSTRUCTION at t
t=1;
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
g=g-mean(s);




%VELOCITY RECONSTRUCTION
z_s=[0:1d-3:20d-2]; %vertical coordinate up to 20cm with 1mm desolution (to be adjusted to match PIV)
u=zeros(length(z_s),length(x_s)-1);w=zeros(length(z_s)-1,length(x_s)-1);
f=fft(s); %Fourrier Modal decomposition
fa=2*abs(f)/length(f); %amplitude of mode
fp=angle(f); %phase of mode
k=[0:(length(s)-1)]/(length(s)-1)*2*pi/dx; %wavenumber of mode
for j=1:length(z_s)
g=0;
h=0;
for i=1:floor(length(s)/2)%/10
g=g-fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*cos(k(i)*x_s(1:end-1)+fp(i))+mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
h=h+fa(i)*sqrt(9.81*k(i))*exp(-z_s(j)*k(i))*sin(k(i)*x_s(1:end-1)+fp(i))-mean(s)/floor(length(s))/2*sqrt(9.81*k(i))*exp(-z_s(j)*k(i));
end
u(j,:)=g;
w(j,:)=h;
end




g=isnan(mean(stack_of_profiles'));
 b=ones(1,ceil(1.7d-2*1/dx))/ceil(1.7d-2*1/dx);
CC=[];
PP=[];

 for i=1:J-1
     if (g(i)+g(i+1))<1 %two in a row without bad profiles
         
         %PHASE SPEED
         s1=delx(i,:); lowpass_s1=filtfilt(b,1,s1); % take a profile then the next one
         s2=delx(i+1,:);lowpass_s2=filtfilt(b,1,s2);
         [c_cov,lags_cov] = xcov(lowpass_s1,lowpass_s2,'coeff');
         c_max = max(c_cov(1:length(s1)));
         lag_off = -lags_cov(find(c_cov==c_max,1,'first'));
         c=dx*lag_off/dt;
         C(i)=c;
         quality_of_phase(i)=c_max;     
         [~,locs] = findpeaks(lowpass_s1);
         a_crest(i)=mean(s1(locs));
         a_crest_error(i)=std(s1(locs));
         
         

%WAVELENGTH & AMPLITUDE
if length(locs)>1
    lambda_crest(i)=mean(diff(locs))*dx;
    lambda_crest_error(i)=std(diff(locs))*dx;
else
    lambda_crest(i)=NaN;
    lambda_crest_error(i)=NaN;
end
[~,locs] = findpeaks(-lowpass_s1);
    a_trough(i)=mean(s1(locs));
    a_trough_error(i)=std(s1(locs));
if length(locs)>1
    lambda_trough(i)=mean(diff(locs))*dx;
    lambda_trough_error(i)=std(diff(locs))*dx;
else
    lambda_trough(i)=NaN;
    lambda_trough_error(i)=NaN;
end
z=find(abs(diff(sign(lowpass_s1)))>1);
if length(z)>1
    lambda_zero(i)=mean(diff(z))*dx*2;
    lambda_zero_error(i)=std(diff(z))*dx*2;
else
    lambda_zero(i)=NaN;
    lambda_zero_error(i)=NaN;
end
       
     else
        C(i)=NaN;
        quality_of_phase(i)=NaN;
        lambda_crest(i)=NaN;
        lambda_crest_error(i)=NaN;
        a_crest(i)=NaN;
        a_crest_error(i)=NaN;
        lambda_trough(i)=NaN;
        lambda_trough_error(i)=NaN;
        a_trough(i)=NaN;
        a_trough_error(i)=NaN;
        lambda_zero(i)=NaN;
        lambda_zero_error(i)=NaN;
         
     end
 end