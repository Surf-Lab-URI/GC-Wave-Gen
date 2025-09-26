function u = medfiltEpsThreshParallele(u,epsilon,thresh,b)
%Westerweel & Scarano (2005): Universal Outlier detection for PIV data
% epsilon=0.1;
% thresh=5;
[J,I]=size(u);
medianres=zeros(J,I);
normfluct=zeros(J,I,2);
% b=2;
eps=0.1;
velcomp=u;
parfor i=1+b:I-b
    for j=1+b:J-b
        neigh=velcomp(j-b:j+b,i-b:i+b);
        neighcol=neigh(:);
        neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
        med=median(neighcol2);
        fluct=velcomp(j,i)-med;
        res=neighcol2-med;
        medianres=median(abs(res));
        normfluct(j,i)=abs(fluct/(medianres+epsilon));
    end
end

info1=(sqrt(normfluct(:,:).^2)>thresh);
u(info1==1)=NaN;
