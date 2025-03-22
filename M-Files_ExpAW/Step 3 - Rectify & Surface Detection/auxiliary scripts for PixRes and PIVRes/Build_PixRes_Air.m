%% Build PixRes_Air

%%% Auxiliary script to build PixRes_Air and make main script lighter

PixRes_Air.XLFV_Surface = XLFV_Surface;
PixRes_Air.LFV_Surface = LFV_Surface;
PixRes_Air.XPIV_LFV_Surface = XPIV_LFV_Surface;
PixRes_Air.PIV_LFV_Surface = PIV_LFV_Surface;
PixRes_Air.PIV_LFV_Surface_smth = filtfilt(ones(1,FiltLength)/FiltLength, 1, PIV_LFV_Surface);
PixRes_Air.PIV_LFV_Surface_smth_phase = angle(hilbert(-PixRes_Air.PIV_LFV_Surface_smth+mean(PixRes_Air.PIV_LFV_Surface_smth,'omitnan')));
PixRes_Air.PIV_Surface = PIV_Surface; % Fits perfectly with PIV1_A

PixRes_Air.pair_index = pair_index;
PixRes_Air.ImageNum_1 = ImageNum_Air1;
PixRes_Air.ImageNum_2 = ImageNum_Air2;
PixRes_Air.PairNum = PairNum;
PixRes_Air.ExpName = ['ExpAW' ExpAW];
PixRes_Air.Acc = Acc;
PixRes_Air.Wind = Wind;
PixRes_Air.Run = runName;

%%%