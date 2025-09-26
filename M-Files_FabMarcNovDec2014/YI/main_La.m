clear
clc

a=['1_03';'1_04';'1_05';'2_01';'2_02';'2_03';'3_01';'3_02';'3_03';'4_01';'4_02';'4_03';'4_05'];


% FOR SURFACE!!!
for n=1:13
A=[];
for i=101:200
    LFV=['G:\FabMarcNovDec2014\Data\Longitudinal\PIVdt10ms_IRlas1_8hz\ExpLCL_',a(n,:),'\PIVRaw\LFV\']
    cd(LFV)
    files=dir(LFV);
    files=files(3:end);
    
    filename=files(i).name;
    load(filename)
    d1=imresize(imgLfv,176.77/103.48);
    imSurf1 = findSurface_simple_ext_force((medfilt2(d1)), 1);
    DX=5.65571e-05;
    t=[1:3499]*DX;
    
    wave_numb=wav_num(t,imSurf1.surface*DX,U10);    % calculate wave number
    A=[A wave_numb];
%     A=[A imSurf1.surface];
end
%     save(['E:\results\wave number\LONG\PIV_lfv\surface',a(n,:),'-101-200.mat'],'A')
end



%WAVE number again 12_22

clear
clc

a=['1_03';'1_04';'1_05';'2_01';'2_02';'2_03';'3_01';'3_02';'3_03';'4_01';'4_02';'4_03';'4_05'];

for n=1:13
	load(['E:\results\wave number\LONG\PIV_lfv\surface\surface',a(n,:),'-101-200.mat'])
    wv_nb=[];
for i=1:100
    surface=A(:,i);
    DX=5.65571e-05;
    t=[1:3499]*DX;   
    wave_numb=wav_num(t,surface);    % calculate wave number
    wv_nb=[wv_nb wave_numb];
end
    save(['E:\results\wave number\LONG\PIV_lfv\wave_number\wave_number',a(n,:),'-101-200.mat'],'wv_nb')
end





% calculate Langmuir number
% load wave rms, that is a.
% load k, that is 
load('E:\results\wave number\LONG\PIV_lfv\exp_wave_num.mat')
load('E:\results\wave\rms\exp1234.mat')
    

rms=[rms1;rms2;rms3;rms4];
t=[101:200]*60/432;
for n=1:4
    M=[];
    for i=1:100        
        temp=rms(n,(i+100)*13-50-13:(i+100)*13-50);
        temp=mean(temp);
        La=Langmuir_number(temp,k3(i));
        M=[M La];
    end
    t1=t+12;
    t2=t+28;
    t3=t+43;
    t4=t+66;
    figure;plot(t3,M)  
    xlabel('time (s)')
    ylabel('Langmuir number')
    title('Langmuir number')
%     saveas(gcf,['E:\results\langmuir number\',num2str(n),'.png'])
end
figure    
subplot(3,1,1)
plot(t4,k4)
   
    ylabel('wave number')
    title('wave number')
    subplot(3,1,2)
    plot(t4,M)  

    ylabel('Langmuir number')
    title('Langmuir number')
    
    subplot(3,1,3)
    DX=5.65571e-05;
    y=[1:511]*DX;
 
    imagesc(t4,y,V_E(:,101:200));
    caxis([0 0.01])
    ylim([0 0.015])
    map=colormap;
    map=sort(map,'descend');
    colormap(map)
     xlabel('time (s)')
     ylabel('depth(m)')
     title('vorticity square')
     saveas(gcf,['E:\results\langmuir number\com_',num2str(n),'.png'])
     
     % VORTICITY
     tt=[1:432]*60/432;
     figure
     imagesc(tt+28,y,V_E);     
%      map=colormap;
%     map=sort(map,'descend');
%     colormap(map)
    colorbar
    xlabel('time(s)')
    ylabel('depth(m)')
    title('vorticity square')
    
    T=[1:635]*IR.DX;
    figure;imagesc(T,T,IR.img)
    colormap gray
    colorbar
    xlabel('length(m) in across wind direction')
    ylabel('length(m) in along wind direction')
    title('IR images when LC just occur')
    saveas(gcf,['E:\results\langmuir number\IR_',num2str(n),'.png'])
    
    
    
    