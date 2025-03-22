% This script loads velocity derivatives from each pair of a selected scene
% and calculates laplacian (first used before function "Lapl.m"

clear
close all
clc

Exp = input('Number of the experiment (put ''4''): ','s');
Scene = input('Number of the scene (put ''3''): ','s');

homePath = 'D:\URI_EXP';
addpath(genpath([homePath,'M_files\']))
resPath = [homePath,'data\EXPW4\EXPW4_Scene3\RESULTS\'];
cartPath = [resPath,'CALCULATED_FIELDS\Cartesian Fields\Cart Fields\'];
cartDiffPath = [resPath,'CALCULATED_FIELDS\Cartesian Fields\CartDiff Fields\'];
transfoPath = [resPath,'Transformations\'];
PairsCartDiff = dir([cartDiffPath,'*cartDiff*']); % find all the pairs in the selected scene

for i = 1:length(PairsCartDiff)
    fname = PairsCartDiff(i).name;
    load([cartDiffPath,fname])
    load([cartPath,'ExpW',Exp,'_Scene',Scene,'_CartVel_',fname(end-7:end)],'PIVRes')
    load([transfoPath,'ExpW',Exp,'_Scene',Scene,'_transfo_',fname(end-7:end)])
    dirX = [0;1]; % directional derivative along x
    dirZ = [1;0]; % directional derivative along z
    L.u_xx = csapsDiff1(cartDiff.u_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of u ONLY in x direction
    L.u_zz = csapsDiff1(cartDiff.u_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of u ONLY in z direction
    L.w_xx = csapsDiff1(cartDiff.w_x, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirX); % Compute 2nd derivative of w ONLY in x direction
    L.w_zz = csapsDiff1(cartDiff.w_z, 0.001, PIVRes.xPIV, PIVRes.zPIV,dirZ); % Compute 2nd derivative of w ONLY in z direction
    L.u_zz = -L.u_zz; % change sign because z is flipped
    L.w_zz = -L.w_zz; % change sign because z is flipped
    % Following lines commented because less efficient (csapsDiff
    % substituted by csapsDiff1)
%     [u_xx2, u_xz2] = csapsDiff(cartDiff.u_x, 0.001, PIVRes.xPIV, PIVRes.zPIV); % Compute 2nd derivative of u in x
%     [u_zx2, u_zz2] = csapsDiff(cartDiff.u_z, 0.001, PIVRes.xPIV, PIVRes.zPIV); % Compute 2nd derivative of u in z
%     [w_xx2, w_xz2] = csapsDiff(cartDiff.w_x, 0.001, PIVRes.xPIV, PIVRes.zPIV); % Compute 2nd derivative of w in x
%     [w_zx2, w_zz2] = csapsDiff(cartDiff.w_z, 0.001, PIVRes.xPIV, PIVRes.zPIV); % Compute 2nd derivative in of w z
%     u_xz2 = -u_xz2; % change sign because z is flipped 
%     u_zz2 = -u_zz2; % change sign because z is flipped
%     w_xz2 = -w_xz2; % change sign because z is flipped
%     w_zz2 = -w_zz2; % change sign because z is flipped
% Note: _xz and _zx must be the same (can be used as a check)
    
    L.Lapl_u = L.u_xx + L.u_zz;
    L.Lapl_w = L.w_xx + L.w_zz;
    L.Lapl_uInt = TransformVelField_decay_0( flipud(L.Lapl_u.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) ); 
    L.Lapl_uInt2 = TransformVelField_decay( flipud(L.Lapl_u.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) ); 
    L.Lapl_wInt = TransformVelField_decay_0( flipud(L.Lapl_w.*cartDiff.Mask), PIVRes, real(transfo.SU(2:end,:)) );
    
    % Use the velocity extrapolated exactly on the surface
    L.u_diff_sum = L.Lapl_uInt(1,:); 
    L.w_diff_sum = L.Lapl_wInt(1,:);
%     % Average of the first two velocity values above the surface on the surface following grid
%     L.u_diff_sum = mean(L.Lapl_uInt(1:2,:)); 
%     L.w_diff_sum = mean(L.Lapl_wInt(1:2,:));
    
    L.DX = CST.DX;
    L.DT = CST.DT;
    
% Following code considering just the first two velocities above the surface (not the first 2 values on the surface following grid as done    
%     Lapl_u = flipud((u_xx+u_zz).*cartDiff.Mask)/CST.DT/CST.DX;% dimensionalise and flip upside down the laplacian field
%     Lapl_w = flipud((w_xx+w_zz).*cartDiff.Mask)/CST.DT/CST.DX;% dimensionalise and flip upside down the laplacian field
%     indLu = zeros(1,size(Lapl_u,2));
%     indLw = zeros(1,size(Lapl_w,2));
%     for ii = 1:size(Lapl_u,2)
%         indLu(ii) = find(~isnan(Lapl_u(:,ii)),1); % find the first node above the surface for Lapl_u
%         indLw(ii) = find(~isnan(Lapl_w(:,ii)),1); % find the first node above the surface for Lapl_w
%     end
%     LindLu = sub2ind([size(Lapl_u,1),size(Lapl_u,2)],indLu,1:size(indLu,2)); % linear indices of the first value above the surface at each x for Lapl_u
%     LindLw = sub2ind([size(Lapl_u,1),size(Lapl_u,2)],indLu,1:size(indLu,2)); % linear indices of the first value above the surface at each x for Lapl_w
%     
%     L.u_diff_sum = mean([Lapl_u(LindLu);Lapl_u(LindLu+1)]); % average of the first two velocities above the surface for Lapl_u
%     L.w_diff_sum = mean([Lapl_w(LindLw);Lapl_w(LindLw+1)]); % average of the first two velocities above the surface for Lapl_w
    
    savePath = [dataPath,'\Lapl Fields\'];
    save([savePath,'ExpW',Exp,'_Scene',Scene,'_Lapl_',Cartesian.PairNum,'_bis.mat'],'L')
    if floor(i/20) == i/20
        disp(num2str(i))
    end
end
