clear
clc

a=['1_03';'1_04';'1_05';'2_01';'2_02';'2_03';'3_01';'3_02';'3_03';'4_01';'4_02';'4_03';'4_05'];
aa=['1_01';'1_02';'1_03';'2_01';'2_02';'2_03';'3_01';'3_02';'3_03';'4_01';'4_02';'4_03'];


for n=1:12
    load(['E:\results\wave number\LONG\PIV_lfv\wave_number\wave_number',a(n,:),'-101-200.mat']) % load wave number
    
    load(['G:\FabMarcNovDec2014\Data\Wg_Pitot\LC',aa(n,:),'\Pitot_WireWG\matFiles\pitot.mat'])   % load pitot tube for wind speed
    load(['E:\results\wave number\LONG\PIV_lfv\surface\surface',a(n,:),'-101-200.mat'])  % to calculate rms of surface
%     load(['E:\results\wave\rms\exp1.mat])   %load rms
% 
    M=[];
    for i=1:100
        k=wv_nb(i);
        if n==1|n==2|n==3;
            s=12;
        elseif n==4|n==5|n==6;
            s=28;
        elseif n==7|n==8|n==9;
            s=43;
        else n==10|n==11|n==12;
            s=56;
        end   
        
        U10=pit_m_s(200*(50+s+(100+i)/432*60)-100:200*(50+s+(100+i)/432*60)+100);
        U10=mean(U10);
        temp=A(:,i);
        DX=5.65571e-05;
        temp=std(temp*DX);
%         temp=rms(n,(i+100)*13-50-13:(i+100)*13-50);
%         temp=mean(temp);                                 
%         to detect the wave gauge to get rms
        La=Langmuir_number(temp,k,U10);
        M=[M La];
    end
    t=[101:140]/432*60;
    
    figure; plot(t,1./M(1:40))
    saveas(gcf,['E:\results\langmuir number\reversed_exp',a(n,:),'.jpg'])
    figure; plot(t,M(1:40))
%     xlabel('wind starting time(s)')
%     ylabel('reverse Langmuir number')
%     title('reverse Langmuir number')
%     save(['E:\results\langmuir number\exp',a(n,:),'.mat'])
    saveas(gcf,['E:\results\langmuir number\exp',a(n,:),'.jpg'])
    save(['E:\results\langmuir number\exp',a(n,:),'.mat'],'M')
end