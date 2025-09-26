clear all
close all


PivRes = 40d-6;
DeltaT= 7e-3;
Sc = PivRes/DeltaT;
yy = [1:100]*8*PivRes;

load('\\Beo\data\ExpLC2_dt7ms_1\ProcessedVelocities\ExpLC2_dt7ms_1_KinEnBudgtTurbFlow.mat')
% load('\\Beo\data\ExpLC2_dt7ms_1\ProcessedVelocities\ExpLC2_dt7ms_1_KinEnBudgtTurbFlow_MoinsOrb.mat')

%% smoothn
Evol = smoothn(KEB.Evol);
Advec = smoothn(KEB.Advec);
TurbTransport = smoothn(KEB.TurbTransport);
MolVisquTransport = smoothn(KEB.MolVisquTransport);
Prod = smoothn(KEB.Prod);
Dissip = smoothn(KEB.Dissip);
Residu = smoothn(KEB.Residu);
KEB.start_image_number

% %% Moyenne glissante
% for i=3:size(KEB.Evol,2)-2
%     Evol(:,i)=(KEB.Evol(:,i-2)+KEB.Evol(:,i-1)+KEB.Evol(:,i)+KEB.Evol(:,i+1)+KEB.Evol(:,i+2))/5;
%     Advec(:,i)=(KEB.Advec(:,i-2)+KEB.Advec(:,i-1)+KEB.Advec(:,i)+KEB.Advec(:,i+1)+KEB.Advec(:,i+2))/5;
%     TurbTransport(:,i)=(KEB.TurbTransport(:,i-2)+KEB.TurbTransport(:,i-1)+KEB.TurbTransport(:,i)+KEB.TurbTransport(:,i+1)+KEB.TurbTransport(:,i+2))/5;
%     MolVisquTransport(:,i)=(KEB.MolVisquTransport(:,i-2)+KEB.MolVisquTransport(:,i-1)+KEB.MolVisquTransport(:,i)+KEB.MolVisquTransport(:,i+1)+KEB.MolVisquTransport(:,i+2))/5;
%     Prod(:,i)=(KEB.Prod(:,i-2)+KEB.Prod(:,i-1)+KEB.Prod(:,i)+KEB.Prod(:,i+1)+KEB.Prod(:,i+2))/5;
%     Dissip(:,i)=(KEB.Dissip(:,i-2)+KEB.Dissip(:,i-1)+KEB.Dissip(:,i)+KEB.Dissip(:,i+1)+KEB.Dissip(:,i+2))/5;
%     Residu(:,i)=(KEB.Residu(:,i-2)+KEB.Residu(:,i-1)+KEB.Residu(:,i)+KEB.Residu(:,i+1)+KEB.Residu(:,i+2))/5;
% end
% 




inter=[1:100];
% n1 = 100;
% n2 = 105;
% n3 = 110;
% n4 = 115;
% n5 = 120;
n1 = 99;
n2 = 104;
n3 = 108;
n4 = 112;
n5 = 119;

n1=104;
n1=111;
n1=115;
n1=117;
%% Evolution
figure, plot(Evol(inter,n1))
hold on, plot(Evol(inter,n2),'c'), plot(Evol(inter,n3),'m'), plot(Evol(inter,n4),'r'), plot(Evol(inter,n5),'k')
%% Advection
figure, plot(Advec(inter,n1))
hold on, plot(Advec(inter,n2),'c'), plot(Advec(inter,n3),'m'), plot(Advec(inter,n4),'r'), plot(Advec(inter,n5),'k')
%% TurbTransport
figure, plot(TurbTransport(inter,n1))
hold on, plot(TurbTransport(inter,n2),'c'), plot(TurbTransport(inter,n3),'m'), plot(TurbTransport(inter,n4),'r'), plot(TurbTransport(inter,n5),'k')
%% MolVisquTransport
figure, plot(MolVisquTransport(inter,n1))
hold on, plot(MolVisquTransport(inter,n2),'c'), plot(MolVisquTransport(inter,n3),'m'), plot(MolVisquTransport(inter,n4),'r'), plot(MolVisquTransport(inter,n5),'k')
%% Prod
figure, plot(Prod(inter,n1))
hold on, plot(Prod(inter,n2),'c'), plot(Prod(inter,n3),'m'), plot(Prod(inter,n4),'r'), plot(Prod(inter,n5),'k')
%% Dissip
figure, plot(Dissip(inter,n1))
hold on, plot(Dissip(inter,n2),'c'), plot(Dissip(inter,n3),'m'), plot(Dissip(inter,n4),'r'), plot(Dissip(inter,n5),'k')
%%Residu
figure, plot(Residu(inter,n1))
hold on, plot(Residu(inter,n2),'c'), plot(Residu(inter,n3),'m'), plot(Residu(inter,n4),'r'), plot(Residu(inter,n5),'k')

%%
n1=115;
figure, plot(Evol(inter,n1)+Advec(inter,n1),-yy)
hold on, plot(TurbTransport(inter,n1),-yy,'c'), plot(MolVisquTransport(inter,n1),-yy,'m'), plot(Prod(inter,n1),-yy,'g'), plot(Dissip(inter,n1),-yy,'k'),  plot(Residu(inter,n1),-yy,'o')

figure, plot(Evol(inter,n1)+Advec(inter,n1),-yy)
figure,plot(TurbTransport(inter,n1),-yy)
figure,plot(MolVisquTransport(inter,n1),-yy)
figure,plot(Prod(inter,n1),-yy)
figure,plot(Dissip(inter,n1),-yy)
figure, plot(Residu(inter,n1),-yy)
