function [FusedPIV] = Pre_process_PIV_Image(FusedPIV);

    clahesize=64;
    numberoftiles1=round(size(FusedPIV,1)/clahesize);
    numberoftiles2=round(size(FusedPIV,2)/clahesize);
    
    FusedPIV=adapthisteq(FusedPIV./max(max(FusedPIV)), 'NumTiles',[numberoftiles1 numberoftiles2], 'ClipLimit', 0.01, 'NBins', 256, 'Range', 'full', 'Distribution', 'uniform');
    FusedPIV=FusedPIV-min(min(FusedPIV));
    FusedPIV=FusedPIV/max(max(FusedPIV))*255;
    
%     h = fspecial('gaussian',32,32);
%     FusedPIV=double(FusedPIV-(imfilter(FusedPIV,h,'replicate')));
%     %FusedPIV=double(FusedPIV-imopen(FusedPIV,strel('disk',16)));
%     FusedPIV=FusedPIV-min(min(FusedPIV));
%     FusedPIV=wiener2(FusedPIV1,[3 3]);
%        
end
