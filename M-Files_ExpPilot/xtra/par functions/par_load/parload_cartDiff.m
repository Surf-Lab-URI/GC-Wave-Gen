function cartDiff = parload_cartDiff(SavePathCartDiff,expName,sceneName,PairNum)

S = load( [SavePathCartDiff 'Movie' expName '_Scene' sceneName '_cartDiff_' PairNum '.mat'] );
cartDiff = S.cartDiff;