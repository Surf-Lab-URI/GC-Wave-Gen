%uniformizer
f=(mean(imgPivsurfcc));
f=f/max(f);
f=1./f;
ff=repmat(f,2048,1);
imgPivsurfcc=imgPivsurfcc.*ff;

%histogram
[h x]=hist(imgPivsurfcc(:),[1:4096]);
c=cumsum(h);
c=c/c(end);

[pks,locs]  = findpeaks(smoothn(h,1000),'MINPEAKHEIGHT',2d4);
m=locs(1)+(locs(1)-find(c>c(locs(1))/10,1,'first'));
M=locs(2)-(find(c>1-(1-c(locs(2)))*0.1,1,'first')-locs(2));
dx=pi/(M-m);
c=[zeros(1,m-1) ((sin(-pi/2:dx:pi/2)+1)/2*4096) ones(1,4096-M)*4096] ;

%find Multiplier
for i=1:2048
    for j=1:2048
M(i,j)=c(round(imgPivsurfcc(i,j)));
    end
end

imgPivsurfcc=imgPivsurfcc.*M;
