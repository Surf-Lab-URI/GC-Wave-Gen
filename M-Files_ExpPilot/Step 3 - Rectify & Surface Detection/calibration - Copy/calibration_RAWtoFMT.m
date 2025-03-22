function calibration_RAWtoFMT(CalDir,SaveFold,Int,MP,FMT)

%% Convert Raw in a FMT format image
for IndIm = Int+1
    eval(['[IM' num2str(IndIm) '] = fliplr(load_Image_IOCoreView_' num2str(MP) 'MP([CalDir(IndIm).folder ''\'' CalDir(IndIm).name]));']);
    ImName = [SaveFold '\' CalDir(IndIm).name(1:end-4) '.' FMT ];
    if strcmp(FMT,'tif')
        FMT2 = 'tiff';
        eval(['imwrite(uint16(IM' num2str(IndIm) '),ImName, ''' FMT2 ''')'])
    elseif strcmp(FMT,'jpg')
        eval(['imwrite(IM' num2str(IndIm) '/512,ImName, ''' FMT ''')'])
    else
        eval(['imwrite(uint16(IM' num2str(IndIm) '),ImName, ''' FMT ''')'])
    end
end
