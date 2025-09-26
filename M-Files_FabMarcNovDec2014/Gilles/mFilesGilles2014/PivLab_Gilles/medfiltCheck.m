function [u,v] = medfiltCheck(u,v,epsilon,thresh)
%Westerweel & Scarano (2005): Universal Outlier detection for PIV data
% epsilon=0.1;
% thresh=5;
[J,I]=size(u);
medianres=zeros(J,I);
normfluct=zeros(J,I,2);
b=1;
eps=0.1;
for c=1:2
    if c==1; velcomp=u;else;velcomp=v;end %#ok<*NOSEM>
    for i=1+b:I-b
        for j=1+b:J-b
            neigh=velcomp(j-b:j+b,i-b:i+b);
            neighcol=neigh(:);
            neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
            med=median(neighcol2);
            fluct=velcomp(j,i)-med;
            res=neighcol2-med;
            medianres=median(abs(res));
            normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
        end
    end
end
info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
u(info1==1)=NaN;
v(info1==1)=NaN;

     