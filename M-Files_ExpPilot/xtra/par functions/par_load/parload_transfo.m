function transfo = parload_transfo(LoadPathTransf,expName,sceneName,PairNum)

S = load( [LoadPathTransf 'Movie' expName '_Scene' sceneName '_transfo_' PairNum '.mat'] );
transfo = S.transfo;