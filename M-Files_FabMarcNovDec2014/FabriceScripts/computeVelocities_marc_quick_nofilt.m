% Function computeVelocities
%%

function compVel =  computeVelocities_marc_quick_nofilt(IM1, IM2, mask1, mask2, IntrWndw, GrdSpc)
% [delx_int dely_int, MASK,cor_mtx,level_mtx]
% tic,
%IM1 = pivIn.IM_a;IM2 = pivIn.IM_b;


IM1=double(IM1).*mask1;
IM2=double(IM2).*mask2;
IM1(isnan(IM1))=nanmean(nanmean(IM1));
IM2(isnan(IM2))=nanmean(nanmean(IM2));
% mask should contain NaN where there's no possible velocity calculation, 1s elsewhere

[h, w] = size(IM1); %image height and width

number_of_levels= length(IntrWndw);

%%
%---------------------- GOBAL CALCULATION ------------------
for lvl = 1:number_of_levels-1; %First level
    
    IW = IntrWndw(lvl); %window size (interrogation window)
    GS = GrdSpc(lvl); %step size (grid spacing)
    
    x = IW/2:GS:(w-IW/2); %x coordinate in image
    y = IW/2:GS:(h-IW/2); %y coordinate in image
    
    bxsNh = floor(1 + (h - IW)/GS); %blocks in height
    bxsNw = floor(1 + (w - IW)/GS); %blocks in width
    
    delx = NaN(bxsNh, bxsNw); % total velocity in x  - accumulates from level to level with INTdel and gdel
    dely = NaN(bxsNh, bxsNw); % total velocity in y  - accumulates from level to level with INTdel and gdel
    dcor = NaN(bxsNh, bxsNw); %
    GLOBAL = NaN(bxsNh, bxsNw); % shows where the GLOBAL Motion is calculated at this level only
    MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
%     IM_MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
    
    
    if lvl==1 %initialize global motion (and correlation) if it's 1st level, otherwise it's passed
        %     from previous level (after being interpolated to the proper grid)
        INTdelx=zeros(bxsNh, bxsNw);
        INTdely=zeros(bxsNh, bxsNw);
        
    end
    
    bxCNTc = 1; %counter initialization (column)
    
    for c = x %loop in column (x, coordinate in image)
        
        bxCNTr = 1; %counter initialization (row)
        
        bdryL = c - IW/2 + 1;  % left and right boundaries of sub-window
        bdryR = c + IW/2;
        
        for r = y %loop in row (y, coordinate in image)
            
            bdryT = r - IW/2 + 1; % top and bottom boundaries of sub-window
            bdryB = r + IW/2;
            
            %global velocity (and correlation) in this block (from previous level)
           
            
            
%           IM_MASK(bxCNTr, bxCNTc)= (~isnan(mask1(r, c))); %if the center of window is not NAN, it's valid
%           VALID_BLOCK = ~isnan(sum(sum(mask(bdryT:bdryB,bdryL:bdryR)))); %if no NANs at all in subwindow
            gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
            gdely_round = round(INTdely(bxCNTr, bxCNTc));
            forwardx=ceil(gdelx_round);
            forwardy=ceil(gdely_round);
            backx=floor(gdelx_round/2)*0;
            backy=floor(gdely_round/2)*0;
            VALID_BLOCK = (~isnan(mask1(r, c))) ; %if center of window not nan
            
            
            
            if VALID_BLOCK %if the PIV can be done (not masked)
                MASK(bxCNTr, bxCNTc)=1;
                try
                    % sub-window blocks
                    
                    bxA = IM1( (bdryT:bdryB) - backy, (bdryL:bdryR) - backx );
                    bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );
                    
                    % sub-window demeaned, windowed, and FFTs
                    bxAmm = bxA-mean(bxA(:));
                    bxBmm = bxB-mean(bxB(:));
                    %                     bxAmm = bxAmm.*Hanning2d;
                    %                     bxBmm = bxBmm.*Hanning2d;
                    fftA = fft2( bxAmm );
                    fftB = fft2( bxBmm );
                    %
                    fftCorr = fftB .* conj(fftA);
                    %PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
                    Xcorr  = fftshift(real( ifft2( fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
                    %[PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in Phase correlation
                    [Xpky, Xpkx]   = find( Xcorr == max(max(Xcorr)) ); %find max in Cross correlation
                    ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                    ldely =  Xpky - IW/2 - 1;
                    T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1)); %3 point gaussian interpolation
                    t = T(:,2);
                    SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                    t = T(2,:);
                    SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                    if isreal(SubpixelY) && isreal(SubpixelX)
                        delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                        dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                        GLOBAL(bxCNTr, bxCNTc) = 1;
                        dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                    else
                        delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                        dely(bxCNTr, bxCNTc) = NaN;
                    end
                    
                end
            
                
            end %end VALID BLOCK CHECK
            bxCNTr = bxCNTr + 1;
            
        end    %end of row for loop
        
        bxCNTc = bxCNTc + 1;
        
    end  %end of column for loop
    %
    % ---------- GLOBAL LEVEL FILTERS -------------------
    %filter & smooth
    %reject displacement more than 1/2 the window size
    %
    %outlier interpolation
    delx=(smoothn(delx,'robust'));
    dely=(smoothn(dely,'robust'));
    %
    %interpolation for next level
    %grid at next level
    IW2 = IntrWndw(lvl+1); %window size (interogation window)
    GS2 = GrdSpc(lvl+1);
    x2 = IW2/2:GS2:(w-IW2/2);
    y2 = IW2/2:GS2:(h-IW2/2);
    %
    %interpolation for next level
    INTdelx = (interp2(x, y', delx, x2, y2', '*spline'));
    INTdely = (interp2(x, y', dely, x2, y2', '*spline'));
end
%
% delx_global=delx;
% dely_global=dely;
% mask_global=MASK;
%
%---------------------- END GLOBAL MOTION ------------------
%%
%------------------- LOCAL (last) LEVEL CALCULATION --------------------------
%
lvl = number_of_levels;
%
IW = IntrWndw(lvl); %window size (interogation window)
GS = GrdSpc(lvl); %step size (grid spacing)
%
x = IW/2:GS:(w-IW/2); %x coordinate in image
y = IW/2:GS:(h-IW/2); %y coordinate in image
%
bxsNh = floor(1 + (h - IW)/GS); %blocks in height
bxsNw = floor(1 + (w - IW)/GS); %blocks in width
%
delx = NaN(bxsNh, bxsNw); %total velocity in x  - accumulates from previous levels with INTdel and gdel
dely = NaN(bxsNh, bxsNw); %total velocity in y  - accumulates from previous levels with INTdel and gdel
dcor = NaN(bxsNh, bxsNw); % total correlation - accumulates from previous levels with INTcor and gcor
MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
CONCa = NaN(bxsNh, bxsNw);
CONCb = NaN(bxsNh, bxsNw);
% IM_MASK = NaN(bxsNh, bxsNw);
%
bxCNTc = 1; %counter initialization
%
for c = x %loop in column (x, coordinate in image)
    
    bxCNTr = 1; %counter initialization (row)
    
    bdryL = c - IW/2 + 1;  % left and right boundaries of sub-window
    bdryR = c + IW/2;
    
    for r = y %loop in row (y, coordinate in image)
        
        bdryT = r - IW/2 + 1; % top and bottom boundaries of sub-window
        bdryB = r + IW/2;
        
%           VALID_BLOCK = ~isnan(sum(sum(mask(bdryT:bdryB,bdryL:bdryR)))); %if no NANs at all in subwindow
            gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
            gdely_round = round(INTdely(bxCNTr, bxCNTc));
            forwardx=ceil(gdelx_round);
            forwardy=ceil(gdely_round);
            backx=floor(gdelx_round/2)*0;
            backy=floor(gdely_round/2)*0;
            VALID_BLOCK = (~isnan(mask1(r, c))) ; %if center of window not nan
        
        if VALID_BLOCK %if the PIV can be done (not masked)
            MASK(bxCNTr, bxCNTc)=1;
            try
                % sub-window blocks
                
                bxA = IM1( (bdryT:bdryB) - backy, (bdryL:bdryR) - backx );
                bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );
                
                  
                % sub-window demeaned, windowed, and FFTs
                bxAmm = bxA-mean(bxA(:));
                bxBmm = bxB-mean(bxB(:));
                
                fftA = fft2( bxAmm );
                fftB = fft2( bxBmm );
                %
                fftCorr = fftB .* conj(fftA);
                %PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
                Xcorr  = fftshift(real( ifft2( fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
                %[PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in Phase correlation
                [Xpky, Xpkx]   = find( Xcorr == max(max(Xcorr)) ); %find max in Cross correlation
                ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                ldely =  Xpky - IW/2 - 1;
                
                T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1)); %3 point gaussian interpolation
                t = T(:,2);
                SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                t = T(2,:);
                SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                if isreal(SubpixelY) && isreal(SubpixelX)
                    delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                    dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                    dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                else
                    delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                    dely(bxCNTr, bxCNTc) = NaN;
                    %delx(bxCNTr, bxCNTc) = gdelx_round + ldelx; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                    %dely(bxCNTr, bxCNTc) = gdely_round + ldely;
                end
            end
            %            
        end %end VALID BLOCK CHECK
        bxCNTr = bxCNTr + 1;
        
    end    %end of row for loop
    
    bxCNTc = bxCNTc + 1;
    
end  %end of column for loop

compVel.INTdelx = INTdelx;
compVel.INTdelz = - INTdely;
compVel.dcor = dcor;
compVel.mask = MASK;
% IM_MASK(IM_MASK==0) = NaN;
% compVel.immask = IM_MASK;
compVel.delta_x = delx; % final result: vectors going downwind (from left to right) are positive
compVel.delta_z = - dely; % final result: vectors going away from surface (upward) are positive

compVel.delx=smoothn(delx,'robust');
compVel.dely=smoothn(-dely,'robust');
%
compVel.xPIV = x; % 0 is upper left corner of image (upwind, away from surface)
compVel.zPIV = y; % 0 is upper left corner of image (upwind, away from surface)
compVel.GS = GrdSpc(end); % final grid spacing
% toc,


