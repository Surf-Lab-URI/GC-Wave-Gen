% test the temperature using all the area in dots. 02/17/2016
clear
clc 
close all
TRAN=['G:\FabMarcNovDec2014\Data\Transverse\PIVdt8ms_IRlas1_8hz\ExpLCTB_1_02\IRMat\'];
files=dir(TRAN)
files=files(3:end,:);
DX=3.870000000000000e-04;  %PIX/METER
DT=0.0231;   % PIX/SECOND


% merge the data of 8 dots
for n=1:8
    load(['E:\results\timeseriesvariance4\old method transverse\data1_02_',num2str(n),'_1-12.mat'])
    loc_max_full_im_x1_1=loc_max_full_im_x1; loc_max_full_im_y1_1=loc_max_full_im_y1;
    load(['E:\results\timeseriesvariance4\old method transverse\data1_02_',num2str(n),'_13-30.mat'])
    loc_max_full_im_x1_2=loc_max_full_im_x1; loc_max_full_im_y1_2=loc_max_full_im_y1;
    load(['E:\results\timeseriesvariance4\old method transverse\data1_02_',num2str(n),'_31-50.mat'])
    loc_max_full_im_x1_3=loc_max_full_im_x1; loc_max_full_im_y1_3=loc_max_full_im_y1;
    loc_x(n,:)=[loc_max_full_im_x1_1 loc_max_full_im_x1_2 loc_max_full_im_x1_3];
    loc_y(n,:)=[loc_max_full_im_y1_1 loc_max_full_im_y1_2 loc_max_full_im_y1_3];
end

old_x=loc_x(:,(24*12+30*15+1):(24*12+30*16));
old_y=loc_y(:,(24*12+30*15+1):(24*12+30*16));
clear loc_max_full_im_x1_1 loc_max_full_im_x1_2 loc_max_full_im_x1_3 loc_max_full_im_x1

close all
tem1=[];tem2=[];
for i=652:650+29
filename=files(i).name
load([TRAN filename])

    figure;imagesc(IR.img)
    colormap gray
    st=std2(IR.img(:,550:600));
    me=mean2(IR.img(:,550:600));
    hold on
    c=contour(medfilt2((IR.img)), me+6*[st,st],'k');
    c=c';
    
    m=find((c(:,1)==me+6*st));
    
    X=zeros(500,500);
    Y=zeros(500,500);
    % X is x, all the areas that you need
    % Y is y, the same as X
    
    for j=1:length(m)-1
        x=c(m(j)+1:m(j+1)-1,1);
        y=c(m(j)+1:m(j+1)-1,2);
        X(j,1:length(x))=x;
        Y(j,1:length(x))=y;   
    end
    
    X(X==0)=NaN; Y(Y==0)=NaN;
    xx=nanmean(X');
    yy=nanmean(Y');
    
    plot(xx,yy,'co')
    xlim([1 635])
    ylim([1 635])
    % next you have got the x and y here which are the area in IR.img that
    % you can compare with the dots location you've had. center of the dots
    % are the location.
    plot(old_y(:,i-651),old_x(:,i-651),'*k')
 
        newpoint = [old_y(:,i-649),old_x(:,i-649)];
        xxyy=[xx;yy]';
        Mdl = KDTreeSearcher(xxyy)
        [n,d] = knnsearch(Mdl,newpoint,'k',1);
        line(xxyy(n,1),xxyy(n,2),'color',[1 0 0],'marker','o','linestyle','none','markersize',10)
    
    x_need=X(n,:);  % we found the area of x -8 dots
    y_need=Y(n,:);  % 8 dots - y in the area
    
    x_n=floor(x_need);
    y_n=floor(y_need);
    plot(x_n,y_n,'*')
   
    xl=1:635;
    yl=1:635;
    [xl,yl]=meshgrid(xl,yl);
    
    for i=1:8
        in=inpolygon(xl,yl,x_n(i,:),y_n(i,:));
        temp=IR.img;
        templ=mean(temp(in==1));
        tem1=[tem1 templ];   % tem1 is the temperature that we need of the mean temperature
    end   
    hold off   
end

tem=[];
for i=1:28   
    tem=[tem tem1(8*(i-1)+1)];
    tem2=[tem tem1(8*(i-1)+2)];
    tem3=[tem tem1(8*(i-1)+3)];
    tem4=[tem tem1(8*(i-1)+4)];
    tem5=[tem tem1(8*(i-1)+5)];
    tem6=[tem tem1(8*(i-1)+6)];
    tem7=[tem tem1(8*(i-1)+7)];
    tem8=[tem tem1(8*(i-1)+8)];
end
a1=(tem(2:end)-tem(1:end-1))/DT;
a2=(tem2(2:end)-tem2(1:end-1))/DT;
a3=(tem3(2:end)-tem3(1:end-1))/DT;
a4=(tem4(2:end)-tem4(1:end-1))/DT;
a5=(tem5(2:end)-tem5(1:end-1))/DT;
a6=(tem6(2:end)-tem6(1:end-1))/DT;
a7=(tem7(2:end)-tem7(1:end-1))/DT;
a8=(tem8(2:end)-tem8(1:end-1))/DT;

figure;plot(a1)
hold on
plot(a2)
plot(a3)
plot(a4)
plot(a5)
plot(a6)
plot(a7)
plot(a8)




