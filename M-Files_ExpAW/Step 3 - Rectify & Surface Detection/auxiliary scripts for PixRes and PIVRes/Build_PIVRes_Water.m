%% Build PIVRes_Water

%%% Auxiliary script to build PIVRes_Water and make main script lighter

PIVRes_Water.BadFramePIVSurfW1 = BadFramePIVSurfW1;
PIVRes_Water.xPIV = CompVelWater.xPIV; % The x coordinates of center of IntrWndws
PIVRes_Water.zPIV = CompVelWater.zPIV; % The y coordinates of center of IntrWndws
PIVRes_Water.GS = CompVelWater.GS; % Final grid spacing
PIVRes_Water.PIVW1_Surface = (PixRes_Water1.PIVW1_Surface(PIVRes_Water.xPIV) )/CompVelWater.GS;%Fits perfectly with imagesc(CompVel.delta_x)
PIVRes_Water.pair_index = pair_index;
PIVRes_Water.ImageNum_1 = ImageNum_Water1;
PIVRes_Water.ImageNum_2 = ImageNum_Water2;
PIVRes_Water.PairNum = PairNum;
PIVRes_Water.ExpName = ['ExpAW' ExpAW];
PIVRes_Water.Acc = Acc;
PIVRes_Water.Wind = Wind;
PIVRes_Water.Run = runName;
PIVRes_Water.PF_Surface = length(PIVRes_Water.zPIV)-PIVRes_Water.PIVW1_Surface+1; % It is needed for transformations;
%it's the surface that would be detected on an upside down PIV image . %Fits perfectly with imagesc(flipud(CompVelWater.delta_x))
I = ismember(PixRes_Water1.XLFV_Water,PIVRes_Water.xPIV);
PIVRes_Water.Phase = PixRes_Water1.LFV_Water_smth_phase(I);

%%%