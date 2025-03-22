function PixRes = parload_PixRes(LoadPathPixRes,expName,sceneName,PairNum)

S = load( [LoadPathPixRes 'Movie' expName '_Scene' sceneName '_Surfaces_' PairNum '.mat'] );
PixRes = S.PixRes;