function CompVel = parload_CompVel(SavePIVPath,expName,sceneName,PairNum)

S = load( [SavePIVPath 'Movie' expName '_Scene' sceneName '_PIV_Velocity_' PairNum '.mat'] );
CompVel = S.CompVel;