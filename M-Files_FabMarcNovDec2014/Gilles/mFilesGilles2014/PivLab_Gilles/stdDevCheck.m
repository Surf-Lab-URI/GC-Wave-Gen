function u = stdDevCheck(u,stdthresh)
% stdthresh=7; 
meanu=nanmean(nanmean(u));
std2u=nanstd(reshape(u,size(u,1)*size(u,2),1));
minvalu=meanu-stdthresh*std2u;
maxvalu=meanu+stdthresh*std2u;
u(u<minvalu)=NaN;
u(u>maxvalu)=NaN;
