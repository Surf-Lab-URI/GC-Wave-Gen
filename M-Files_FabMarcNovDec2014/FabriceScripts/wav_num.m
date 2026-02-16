function wave_numb=wav_num(t,x1,DX)


sp=spectrum(x1,length(x1),0,hanning(length(x1)));
sp=sp(:,1);
k=[0:length(sp)-1]/(length(sp)-1)*(2*pi/2/DX);
[y_max index]=max(sp);
wave_numb=k(index);
end
