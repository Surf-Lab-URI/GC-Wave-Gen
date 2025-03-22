function [PIVSurfA_Undistorted] = PIVSurfA_LFV_LensDistCorr(IM1,cameraParams)

PIVSurfA_Undistorted = undistortImage(IM1,cameraParams);
