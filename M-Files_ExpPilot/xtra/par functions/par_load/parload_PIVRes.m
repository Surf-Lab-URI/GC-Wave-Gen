function PIVRes = parload_PIVRes(LoadPathPixRes,expName,sceneName,PairNum)

S = load( [LoadPathPixRes 'Movie' expName '_Scene' sceneName '_Surfaces_' PairNum '.mat'] );
PIVRes = S.PIVRes;