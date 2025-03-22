function [FusedPIV] = Pre_process_PIV_Image_Water_IM1(FusedPIV1)

    clahesize=64;
    numberoftiles1=round(size(FusedPIV1,1)/clahesize);
    numberoftiles2=round(size(FusedPIV1,2)/clahesize);
    % ClipLimit = 0.02;%0.002;
    ClipLimit = 0.02;
    
    FusedPIV=adapthisteq(FusedPIV1./max(max(FusedPIV1)), 'NumTiles',[numberoftiles1 numberoftiles2], 'ClipLimit', ClipLimit, 'NBins', 1028, 'Range', 'full', 'Distribution', 'uniform');
    FusedPIV=FusedPIV-min(min(FusedPIV));
    FusedPIV=FusedPIV/max(max(FusedPIV))*255;
    
%     h = fspecial('gaussian',32,32);
%     FusedPIV=double(FusedPIV-(imfilter(FusedPIV,h,'replicate')));
%     %FusedPIV=double(FusedPIV-imopen(FusedPIV,strel('disk',16)));
%     FusedPIV=FusedPIV-min(min(FusedPIV));
%     FusedPIV=wiener2(FusedPIV1,[3 3]);
%        
end
