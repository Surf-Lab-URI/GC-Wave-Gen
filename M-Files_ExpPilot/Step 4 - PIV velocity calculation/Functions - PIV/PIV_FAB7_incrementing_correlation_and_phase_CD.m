function [CompVel] = PIV_FAB7_incrementing_correlation_and_phase_CD(IM1, IM2, Mask, IntrWndw, GrdSpc)


tic
IM1 = double(IM1) .* Mask;
IM2 = double(IM2) .* Mask;
IM1(isnan(IM1)) = nanmean(nanmean(IM1));
IM2(isnan(IM2)) = nanmean(nanmean(IM2));
IM1_D=IM1;
IM2_D=IM2;
% mask should contain NaN where there's no possible velocity calculation, 1s elsewhere

[h, w] = size(IM1); %image height and width
[X1,Y1] = meshgrid([1:w], [1:h]);

% IntrWndw=[64 32 16 8];
% GrdSpc=[32 16 8 4];
number_of_levels= length(IntrWndw);
% The first four, five, six levels are called the "Global Level" and the last level is called the "Local Level".

%%
%---------------------- GOBAL CALCULATION ------------------

for lvl = 1:number_of_levels-1; %First level
      
IW = IntrWndw(lvl); %window size (interogation window)
GS = GrdSpc(lvl); %step size (grid spacing)

x = IW/2:GS:(w-IW/2); %x coordinate in image
y = IW/2:GS:(h-IW/2); %y coordinate in image

bxsNh = floor(1 + (h - IW)/GS); %blocks in height
bxsNw = floor(1 + (w - IW)/GS); %blocks in width

delx = NaN(bxsNh, bxsNw); % total velocity in x  - accumulates from level to level with INTdel and gdel
dely = NaN(bxsNh, bxsNw); % total velocity in y  - accumulates from level to level with INTdel and gdel
dcor = NaN(bxsNh, bxsNw); % total correlation  - accumulates from level to level with INTcor and gcor
GLOBAL = NaN(bxsNh, bxsNw); % shows where the GLOBAL Motion is calculated at this level only
MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
dwhatlevel = NaN(bxsNh, bxsNw);

if lvl==1 %initialize global motion (and correlation) if it's 1st level, otherwise it's passed from previous level (after beeing interpolated to the proper grid)
INTdelx=zeros(bxsNh, bxsNw); 
INTdely=zeros(bxsNh, bxsNw);  
INTdcor=zeros(bxsNh, bxsNw); 
INTdwhatlevel=zeros(bxsNh, bxsNw);  
end

bxCNTc = 1; %counter initialization (column)

for c = x %loop in column (x)
    
    bxCNTr = 1; %counter initialization (row)
    
    bdryL = c - IW/2 + 1;  % left and right boundaries of sub-window
    bdryR = c + IW/2;
    
    for r = y %loop in row (y)
        
        bdryT = r - IW/2 + 1; % top and bottom boundaries of sub-window
        bdryB = r + IW/2;
        
        %global velocity (and correlation) in this block (from previous level)
        gdelx = INTdelx(bxCNTr, bxCNTc); 
        gdely = INTdely(bxCNTr, bxCNTc);
        gcor = INTdcor(bxCNTr, bxCNTc);
        gwhatlevel = INTdwhatlevel(bxCNTr, bxCNTc);
        VALID_BLOCK=(~isnan(Mask(r, c))); %if the center of window is not NAN, it's valid
       
            
                       
        if VALID_BLOCK %if the PIV can be done (not masked)
           MASK(bxCNTr, bxCNTc)=1;
       
            % sub-window blocks
            bxA = IM1_D( (bdryT:bdryB), (bdryL:bdryR));
            bxB = IM2_D( (bdryT:bdryB), (bdryL:bdryR));
                       
            %hanning sub-window
            N=size(bxA);
            Nx=N(2);  
            Ny=N(1);
            Hanning2d=0.5*(1-cos(2*pi*[0:Ny-1]'/(Ny-1)))*0.5*(1-cos(2*pi*[0:Nx-1]/(Nx-1)));

            % sub-window demeaned, windowed, and FFTs
            bxAmm = bxA-mean(bxA(:));
            bxBmm = bxB-mean(bxB(:));
            bxAmm = bxAmm.*Hanning2d;
            bxBmm = bxBmm.*Hanning2d;
            fftA = fft2( bxAmm );
            fftB = fft2( bxBmm );
            fftCorr = fftB .* conj(fftA);
            PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
            Xcorr=fftshift(real(ifft2(fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
            [PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in Phase correlation
            [Xpky, Xpkx] = find( Xcorr == max(max(Xcorr)) ); %find max in Cross correlation
             
              if ( length(PHpky)==1 && length(PHpkx)==1 && length(Xpky)==1 && length(Xpkx)==1 && abs(Xpky-PHpky)<2 && abs(Xpkx-PHpkx)<2  ) %if local calculation OK. (both phase and cross correlation agree with single peak)
                ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                ldely = Xpky - IW/2 - 1;
                dcor(bxCNTr, bxCNTc) = Xcorr(Xpky, Xpkx); %local correlation calculated here
                %if (dcor(bxCNTr, bxCNTc)>INTdcor(bxCNTr, bxCNTc)-0.5/lvl && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                if (dcor(bxCNTr, bxCNTc)>0.5 && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                    GLOBAL(bxCNTr, bxCNTc)=1; %a Global velocity estimate was made at this level
                    dwhatlevel(bxCNTr, bxCNTc) = lvl; 
                
                   [Xcorr_y Xcorr_x]=size(Xcorr);
                   if ((Xpkx <= Xcorr_x-1) && (Xpky <= Xcorr_y-1) && (Xpkx >= 2) && (Xpky >= 2)) %is subpix interpolation possible (enough room around peak?)
                   T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1)); %3 point gaussian interpolation
                   t = T(:,2);
                   SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                   t = T(2,:);
                   SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                                                         
                       if (isreal([SubpixelY SubpixelX]) && abs(SubpixelY)<1 && abs(SubpixelX)<1) %valid subpixel
                       delx(bxCNTr, bxCNTc) = gdelx + ldelx + SubpixelX;; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                       dely(bxCNTr, bxCNTc) = gdely + ldely + SubpixelY;; % another possibility is delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + (SubpixelX + (gdelx-gdelx_round))/2 ;
                       else % subpixel not good -> keep just local velocity
                       delx(bxCNTr, bxCNTc) = gdelx + ldelx; %gdelx potentially has a subpix estimate from previous level
                       dely(bxCNTr, bxCNTc) = gdely + ldely;
                       end
                    
                   else % subpix estimate not possible -> keep just local velocity
                   delx(bxCNTr, bxCNTc) = gdelx + ldelx; %gdelx potentially has a subpix estimate from previous level
                   dely(bxCNTr, bxCNTc) = gdely + ldely;
                   end %end subpix possible
                else %if correlation NOT better than before, keep old
                delx(bxCNTr, bxCNTc) = gdelx;
                dely(bxCNTr, bxCNTc) = gdely; 
                dcor(bxCNTr, bxCNTc) = gcor;
                dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel;
                end
                                                  
                else %no local calculation possibble -> keep global velocity and correlation from previous level
                %total velocity - global
                delx(bxCNTr, bxCNTc) = gdelx;  %gdelx potentially has a subpix estimate from previous level
                dely(bxCNTr, bxCNTc) = gdely; 
                dcor(bxCNTr, bxCNTc) = gcor;
                dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel;
             end %local calculation possibble
                
            
        end %end VALID BLOCK CHECK
        bxCNTr = bxCNTr + 1;
        
    end    %end of row for loop
    
    bxCNTc = bxCNTc + 1;
    
end  %end of column for loop

% ---------- GLOBAL LEVEL FILTERS -------------------
%filter & smooth 
%reject displacement more than 1/2 the window size
 if lvl==1 %if level=1, do not keep the 0s from initialization where the calculation failed! instead replace with NaNs
 delx=delx.*GLOBAL;
 dely=dely.*GLOBAL;
 dcor=dcor.*GLOBAL;
 dwhatlevel=dwhatlevel.*GLOBAL;
 end
 
    delx = (smoothn(delx,'robust'));
    dely = (smoothn(dely,'robust'));
       
    [X,Y] = meshgrid(x, y);
    %Interpolation of velocity to image resolution
    U1 = interp2(X,Y,delx,X1,Y1,'*spline'); 
    V1 = interp2(X,Y,dely,X1,Y1,'*spline');
    
    %Warping both images according to velocity (centered difference)
    IM1_D= interp2(1:size(IM1,2),(1:size(IM1,1))',IM1,X1-U1/2,Y1-V1/2,'*linear');
    IM2_D= interp2(1:size(IM2,2),(1:size(IM2,1))',IM2,X1+U1/2,Y1+V1/2,'*linear');


    %interpolation for next level
    %grid at next level
    IW2 = IntrWndw(lvl+1); %window size (interogation window)
    GS2 = GrdSpc(lvl+1);
    x2 = IW2/2:GS2:(w-IW2/2);
    y2 = IW2/2:GS2:(h-IW2/2);
    %
    %interpolation for next level
    
INTdelx = U1(y2,x2);
INTdely = V1(y2,x2);
INTdcor = (interp2(x, y, dcor, x2, y2', '*nearest')); %neareast better?
%INTdcor(isnan(INTdcor))=0;
INTdwhatlevel = (interp2(x, y, dwhatlevel, x2, y2', '*nearest')); %neareast better?
%INTdwhatlevel(isnan(INTdwhatlevel))=0;

end


%---------------------- END GLOBAL MOTION ------------------
%%
%------------------- LOCAL (last) LEVEL CALCULATION --------------------------

lvl = number_of_levels;
    
IW = IntrWndw(lvl); %window size (interogation window)
GS = GrdSpc(lvl); %step size (grid spacing)

x = IW/2:GS:(w-IW/2); %x coordinate in image
y = IW/2:GS:(h-IW/2); %y coordinate in image

bxsNh = floor(1 + (h - IW)/GS); %blocks in height
bxsNw = floor(1 + (w - IW)/GS); %blocks in width

delx = NaN(bxsNh, bxsNw); %total velocity in x  - accumulates from previous levels with INTdel and gdel
dely = NaN(bxsNh, bxsNw); %total velocity in y  - accumulates from previous levels with INTdel and gdel
dcor = NaN(bxsNh, bxsNw); % total correlation - accumulates from previous levels with INTcor and gcor
sub = NaN(bxsNh, bxsNw); % where subpixel correlation was used at this level
LOCAL = NaN(bxsNh, bxsNw); % where delx and dely values were modified at this level (others are leftover passthrough from previous levels)
MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
corrmtx = NaN(bxsNh, bxsNw); %correlation matrix  - at this level only
dwhatlevel = NaN(bxsNh, bxsNw); 

bxCNTc = 1; %counter initialization

for c = x %loop in column (x)
    
    bxCNTr = 1; %counter initialization
    
    bdryL = c - IW/2 + 1;  % left and right boundaries of sub-window
    bdryR = c + IW/2;
    
    for r = y %loop in row (y)
        
        bdryT = r - IW/2 + 1; % top and bottom boundaries of sub-window
        bdryB = r + IW/2;
        
        %global velocity (and correlation) in this block (from previous level)
        gdelx = INTdelx(bxCNTr, bxCNTc); 
        gdely = INTdely(bxCNTr, bxCNTc);
        gcor = INTdcor(bxCNTr, bxCNTc);
        gwhatlevel = INTdwhatlevel(bxCNTr, bxCNTc);
        VALID_BLOCK=(~isnan(Mask(r, c))); %if the center of window is not NAN, it's valid
               
        if VALID_BLOCK %if the PIV can be done (not masked)
            MASK(bxCNTr, bxCNTc)=1;
           
                
            % sub-window blocks
            bxA = IM1_D( (bdryT:bdryB), (bdryL:bdryR));
            bxB = IM2_D( (bdryT:bdryB), (bdryL:bdryR));
            
            %hanning sub-window 
            N=size(bxA);
            Nx=N(2);  
            Ny=N(1);
            Hanning2d=0.5*(1-cos(2*pi*[0:Ny-1]'/(Ny-1)))*0.5*(1-cos(2*pi*[0:Nx-1]/(Nx-1)));

            % sub-window  demeaned, FFTs
            bxAmm = bxA-mean(bxA(:));
            bxBmm = bxB-mean(bxB(:));
            bxAmm = bxAmm.*Hanning2d;
            bxBmm = bxBmm.*Hanning2d;
            fftA = fft2( bxAmm );
            fftB = fft2( bxBmm );
            fftCorr = fftB .* conj(fftA);
            PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
            Xcorr=fftshift(real(ifft2(fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
            [PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in phase correlation
            [Xpky, Xpkx] = find( Xcorr == max(max(Xcorr)) );%find max in cross correlation
            
               if ( length(PHpky)==1 && length(PHpkx)==1 && length(Xpky)==1 && length(Xpkx)==1 && Xpky==PHpky && Xpkx==PHpkx  ) %if local calculation OK. (both phase and cross correlation agree with single peak)
                ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                ldely = Xpky - IW/2 - 1;
                dcor(bxCNTr, bxCNTc) = Xcorr(Xpky, Xpkx); %local correlation calculated here
                corrmtx(bxCNTr, bxCNTc)= Xcorr(Xpky, Xpkx);
                if (dcor(bxCNTr, bxCNTc)>INTdcor(bxCNTr, bxCNTc)-0.5/lvl && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                %if (dcor(bxCNTr, bxCNTc)>0 && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before

                LOCAL(bxCNTr, bxCNTc)=1;
                dwhatlevel(bxCNTr, bxCNTc) = lvl;  
                
                   [Xcorr_y Xcorr_x]=size(Xcorr);
                   if ((Xpkx <= Xcorr_x-1) && (Xpky <= Xcorr_y-1) && (Xpkx >= 2) && (Xpky >= 2)) %is subpix interpolation possible
                   T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1));
                   t = T(:,2);
                   SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                   t = T(2,:);
                   SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                                                         
                       if (isreal([SubpixelY SubpixelX]) && abs(SubpixelY)<1 && abs(SubpixelX)<1) %valid subpixel is good
                       delx(bxCNTr, bxCNTc) = gdelx + ldelx + SubpixelX;
                       dely(bxCNTr, bxCNTc) = gdely + ldely + SubpixelY;
                       sub(bxCNTr, bxCNTc) = 1;
                       else % subpixel not good -> keep just local velocity
                       delx(bxCNTr, bxCNTc) = gdelx + ldelx;
                       dely(bxCNTr, bxCNTc) = gdely + ldely;
                       end
                    
                   else % subpix estimate not possible -> keep just local velocity
                   delx(bxCNTr, bxCNTc) = gdelx + ldelx;
                   dely(bxCNTr, bxCNTc) = gdely + ldely;
                   end %end subpix possible
                else %if correlation NOT better than before, keep old
                delx(bxCNTr, bxCNTc) = gdelx;
                dely(bxCNTr, bxCNTc) = gdely; 
                dcor(bxCNTr, bxCNTc) = gcor;
                dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel; 
                end %end correlation NOT better than before
                                                                  
                else %no local calculation possible -> keep global velocity
                %total velocity - global
                delx(bxCNTr, bxCNTc) = gdelx;
                dely(bxCNTr, bxCNTc) = gdely; 
                dcor(bxCNTr, bxCNTc) = gcor;
                dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel; 
                end %local calculation possibble
                
            
        
        end %end VALID BLOCK CHECK
        bxCNTr = bxCNTr + 1;
        
    end    %end of row for loop
    
    bxCNTc = bxCNTc + 1;
end

    
% --------------------- END LOCAL LEVEL CALCULATION --------------------- %
%%


CompVel.INTdelx = INTdelx;
CompVel.INTdelz = - INTdely;
CompVel.dcor = dcor;
CompVel.Mask = MASK;

CompVel.dwhatlevel=dwhatlevel;
CompVel.last_level_corr=corrmtx;
CompVel.subpixel=sub;
CompVel.last_level_vel=LOCAL;

CompVel.delta_x = delx; % Final Result: vectors going downwind (from left to right) are positive
CompVel.delta_z = - dely; % Final Result: vectors going away from surface (upward) are positive

CompVel.xPIV = x; % 0 is upper left corner of image (upwind, away from surface)
CompVel.zPIV = y; % 0 is upper left corner of image (upwind, away from surface)
CompVel.GS = GrdSpc(end); % Final grid spacing
CompVel.IW = IntrWndw(end);
toc


