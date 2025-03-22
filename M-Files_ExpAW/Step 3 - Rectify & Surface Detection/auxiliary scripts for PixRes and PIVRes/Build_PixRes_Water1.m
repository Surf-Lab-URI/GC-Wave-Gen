%% Build PixRes_Water1

%%% Auxiliary script to build PixRes_Water1 and make main script lighter

PixRes_Water1.BadFramePIVSurfW1 = BadFramePIVSurfW1;
PixRes_Water1.XPIVW_PIVSurfW1_Surface = XPIVW_PIVSurfW1_Surface; % PIVWater x-axis in PIV Water coordinates
PixRes_Water1.PIVW_PIVSurfW1_Surface = PIVW_PIVSurfW1_Surface; % PIVWater surface in PIV Water coordinates
PixRes_Water1.PIVW_PIVSurfW1_Surface_smth = filtfilt(ones(1,round(FiltLength/3.333))/(round(FiltLength/3.333)), 1, PIVW_PIVSurfW1_Surface); % smoothed version to calculate gradients
PixRes_Water1.PIVW1_Surface = PIVW1_Surface; % Fits perfectly with imagesc(PIV1_W)

[XPIVW_LFV_Surface,PIVW_LFV_Surface_smth] = transform_phase_from_PIVAir_to_PIVWater(XPIV_LFV_Surface,PixRes_Air.PIV_LFV_Surface_smth); %%% Find phase from LFV in PIVWater coordinates
% This is needed to calculate the phase in water space
% coordinates with long components from LFV

PixRes_Water1.XLFV_Water = XPIVW_LFV_Surface;
PixRes_Water1.LFV_Water_smth = PIVW_LFV_Surface_smth;
PixRes_Water1.LFV_Water_smth_phase = angle(hilbert(-PixRes_Water1.LFV_Water_smth+mean(PixRes_Water1.LFV_Water_smth,'omitnan')));

PixRes_Water1.pair_index = pair_index;
PixRes_Water1.ImageNum_1 = ImageNum_Water1;
PixRes_Water1.PairNum = PairNum;
PixRes_Water1.ExpName = ['ExpAW' ExpAW];
PixRes_Water1.Acc = Acc;
PixRes_Water1.Wind = Wind;
PixRes_Water1.Run = runName;

%%%