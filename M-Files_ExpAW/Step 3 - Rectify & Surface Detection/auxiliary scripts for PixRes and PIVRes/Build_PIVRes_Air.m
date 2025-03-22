%% Build PIVRes_Air

%%% Auxiliary script to build PixRes_Air and make main script lighter

PIVRes_Air.PIV_Surface_PIVAir = PIV_Surface_PIVAir; % this is the surface matching CompVelWater length (that is a bit larger than CompVelAir) in PIVAir coordinates
PIVRes_Air.BadFramePIVSurfLFV = BadFramePIVSurfLFV;
PIVRes_Air.xPIV = CompVelAir.xPIV; % The x coordinates of center of IntrWndws
PIVRes_Air.zPIV = CompVelAir.zPIV; % The y coordinates of center of IntrWndws
PIVRes_Air.GS = CompVelAir.GS; % Final grid spacing
PIVRes_Air.PIV_Surface = (PixRes_Air.PIV_Surface(PIVRes_Air.xPIV) )/CompVelAir.GS;%Fits perfectly with imagesc(CompVel.delta_x)
PIVRes_Air.pair_index = pair_index;
PIVRes_Air.ImageNum_1 = ImageNum_Air1;
PIVRes_Air.ImageNum_2 = ImageNum_Air2;
PIVRes_Air.PairNum = PairNum;
PIVRes_Air.ExpName = ['ExpAW' ExpAW];
PIVRes_Air.Acc = Acc;
PIVRes_Air.Wind = Wind;
PIVRes_Air.Run = runName;
PIVRes_Air.PF_Surface = length(PIVRes_Air.zPIV)-PIVRes_Air.PIV_Surface+1; % It is needed for transformations;
%it's the surface that would be detected on an upside down PIV image . %Fits perfectly with imagesc(flipud(CompVelAir.delta_x))
I = ismember(PixRes_Air.XPIV_LFV_Surface,PIVRes_Air.xPIV);
PIVRes_Air.Phase = PixRes_Air.PIV_LFV_Surface_smth_phase(I);

%%%