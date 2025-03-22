%% Pitot
clear
% close all
clc

% Find files
LoadPath = '\\spray3\d\data\Pitot\Pitot_tests\';
Acc = ''; %'acc0.07';
WS = 'WS7';
Dir = dir([LoadPath '04052023*' Acc '_*' WS '_*.dat']);
dt = 0.02; % dt for datarate=50Hz
%%% Parameters for 03/27/2023
% Ramp = [20 40 80];
% WindSpeed = [4:8];
%%% Parameters for 04/05/2023
% Ramp = [ ... ];
% WindSpeed = [5:8];

for i = 1:length(Dir) %1:length(Dir)
    
    % Load Pitot File
    PitotName = Dir(i).name;
    fclose('all');
    fid = fopen([LoadPath PitotName]);
    headerlength = 1;
    header = fread(fid, headerlength, 'double', 'ieee-be');
    PitotDat = fread(fid,inf,'double', 'ieee-be');
    fclose(fid);
    
    % Pitot Downstream and Upstream in Volts
    PitotDOWN_V = PitotDat(1:length(PitotDat)/2);
    PitotUP_V = PitotDat(length(PitotDat)/2+1:end);
    
    %%% Conversion to m/s (from \\spray3\d\data\Pitot\Calibration_Pitot.xlsx)
    % Downstream
    WCdown_inch = 0.1994*PitotDOWN_V-0.0081;
    Pdown_Pa = WCdown_inch*9806*0.0254;
    VelDown{i} = sqrt(2*Pdown_Pa/1.2); % rho_air = 1.2 kg/m3
    %Upstream
    WCup_inch = 0.1997*PitotUP_V-0.0078;
    Pup_Pa = WCup_inch*9806*0.0254;
    VelUp{i} = sqrt(2*Pup_Pa/1.2); % rho_air = 1.2 kg/m3
    
end

%%% Uncomment this only for tests of 03/27/2023 
% VelUp(:,16) = [];
% VelDown(:,16) = [];

%%
Ncond = 15; % for tests of 03/27/2023
Ncond = length(Dir); % for tests of 04/05/2023
figure;hold on
ylabel('velocity (m/s)');
xlabel('time (s)')
title('Upstream Pitot Tube')
for i = 1:Ncond
    count = 0;
    VUp = [];
    VDown = [];
    
    %%% Use this for tests of 03/27/2023
    % for ii = (i-1)*5+1:i*5
    %     count = count+1;
    %     VUp(:,count) = VelUp{ii};
    %     VDown(:,count) = VelDown{ii};
    % end
    % velUP{i} = mean(VUp,2);
    % velDOWN{i} = mean(VDown,2);
    
    velUP{i} = VelUp{i};
    velDOWN{i} = VelDown{i};
    
    % Plot results upstream
    t = (0:length(velUP{i})-1)*dt;
    plot(t(1:100:end),velUP{i}(1:100:end),'r')
end

% figure;hold on
% ylabel('velocity (m/s)');
% xlabel('time (s)')
% title('Downstream Pitot Tube')
% for i = 1:Ncond
%     % Plot results downstream
%     t = (0:length(velDOWN{i})-1)*dt;
%     plot(t,velDOWN{i},'b')
% end