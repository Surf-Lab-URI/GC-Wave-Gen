% Function computeVelocities
%% Object:
% computes velocities in the air above the water (ones in mask), from
% particle displacement in fused images
%% Arguments: 
% experiment name, image pair number, fused image a, fused image b, mask of
% ones and zeros
%% Result: 
% experiment name, image pair number, delta_x (horizontal displacements),
% delta_z (vertical displacements), correlations, last level reached, mask
% used
%% Author:
% Fabrice Veron, modified by Marc Buckley
%% Last update:
% 05/08/2013
%%
function compVel =  computeVelocities_LC1(IM1, IM2, mask)
% [delx_int dely_int, MASK,cor_mtx,level_mtx]
IM1=double(IM1);
IM2=double(IM2);
% mask should contain NaN where there's no possible velocity calculation, 1s elsewhere

[h, w] = size(IM1); %image height and width

IntrWndw=[128 64 32 16];
GrdSpc=[64 32 16 8];
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
    dcor = NaN(bxsNh, bxsNw); % total correlation  - accumulates from level to level with INTcor and gcor
    GLOBAL = NaN(bxsNh, bxsNw); % shows where the GLOBAL Motion is calculated at this level only
    MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
    dwhatlevel = NaN(bxsNh, bxsNw);

    if lvl==1 %initialize global motion (and correlation) if it's 1st level, otherwise it's passed
        %     from previous level (after being interpolated to the proper grid)
        INTdelx=zeros(bxsNh, bxsNw);
        INTdely=zeros(bxsNh, bxsNw);
        INTdcor=zeros(bxsNh, bxsNw);
        INTdwhatlevel=zeros(bxsNh, bxsNw);
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
            gdelx = INTdelx(bxCNTr, bxCNTc);
            gdely = INTdely(bxCNTr, bxCNTc);
            gcor = INTdcor(bxCNTr, bxCNTc);
            gwhatlevel = INTdwhatlevel(bxCNTr, bxCNTc);
            VALID_BLOCK = (~isnan(mask(r, c))); %if the center of window is not NAN, it's valid
            gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
            gdely_round = round(INTdely(bxCNTr, bxCNTc));
            
            if VALID_BLOCK %if the PIV can be done (not masked)
                MASK(bxCNTr, bxCNTc)=1;
                if ( (bdryT + gdely_round )>=1 && (bdryB + gdely_round)<=h && ...
                        (bdryL + gdelx_round )>=1 && (bdryR + gdelx_round)<=w ) 
                % check that flow points within the boundary of 2nd image
                    
                    % sub-window blocks
                    bxA = IM1(bdryT:bdryB, bdryL:bdryR);
                    bxB = IM2( (bdryT:bdryB) + gdely_round, (bdryL:bdryR) + gdelx_round );
                    
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
                    %
                    fftCorr = fftB .* conj(fftA);
                    PHcorr = fftshift( abs( ifft2( fftCorr./(abs(fftCorr)+eps) ) ) ); %Phase correlation
                    Xcorr  = fftshift(real( ifft2( fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
                    [PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in Phase correlation
                    [Xpky, Xpkx]   = find( Xcorr == max(max(Xcorr)) ); %find max in Cross correlation
                    %
                    if ( length(PHpky)==1 && length(PHpkx)==1 && length(Xpky)==1 && length(Xpkx)==1 && Xpky==PHpky && Xpkx==PHpkx  ) 
                    %if local calculation OK. (both maxs of phase and cross correlation are at the same location (same peak))
                        ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                        ldely =  Xpky - IW/2 - 1;
                        dcor(bxCNTr, bxCNTc) = Xcorr(Xpky, Xpkx); %local correlation calculated here
                        if (dcor(bxCNTr, bxCNTc)>0 && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                        %if (dcor(bxCNTr, bxCNTc)>INTdcor(bxCNTr, bxCNTc)-0.5/lvl && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                        %if (dcor(bxCNTr, bxCNTc)>0.75 && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                            GLOBAL(bxCNTr, bxCNTc) = 1; %a Global velocity estimate was made at this level
                            dwhatlevel(bxCNTr, bxCNTc) = lvl;
                            [Xcorr_y Xcorr_x] = size(Xcorr);
                            if ((Xpkx <= Xcorr_x-1) && (Xpky <= Xcorr_y-1) && (Xpkx >= 2) && (Xpky >= 2)) 
                            %is subpix interpolation possible (enough room around peak?)
                                T=log(Xcorr(Xpky-1:Xpky+1,Xpkx-1:Xpkx+1)); %3 point gaussian interpolation
                                t = T(:,2);
                                SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                                t = T(2,:);
                                SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                                                         
                                if (isreal([SubpixelY SubpixelX]) && abs(SubpixelY)<1 && abs(SubpixelX)<1) %valid subpixel
                                    delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                    dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY; % another possibility is delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + (SubpixelX + (gdelx-gdelx_round))/2 ;
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
                    
                    
                else %if current boxes aren't within image boundaries - no velocity measured at level
                    %total velocity - global
                    delx(bxCNTr, bxCNTc) = gdelx;
                    dely(bxCNTr, bxCNTc) = gdely;
                    dcor(bxCNTr, bxCNTc) = gcor;
                    dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel;
                end %end of boundary check if statement
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

    %Local filter - reject displacement locally where the fluctuation is more than twice the local std
    b=1; %3x3 filter; b=2 -> 5x5 filter etc...
    delx=padarray(delx,[b b],'replicate','both');
    dely=padarray(dely,[b b],'replicate','both');
    valid=~isnan(padarray(MASK,[b b],'replicate','both'));
    [J,I]=size(delx);
    epsilon=0.1*(1/lvl); %becomes more and more stringent with level.
    thresh=2;
    median_Unbrs=NaN(J,I);
    median_Vnbrs=NaN(J,I);
    median_u=NaN(J,I);
    median_v=NaN(J,I);
    fluct_u=NaN(J,I);
    fluct_v=NaN(J,I);
    normfluct_u=NaN(J,I);
    normfluct_v=NaN(J,I);
    for i=1+b:I-b
        for j=1+b:J-b
            if valid(j,i)
                u=delx(j-b:j+b,i-b:i+b);
                v=dely(j-b:j+b,i-b:i+b);
                neighcol_u=u(:);
                neighcol_v=v(:);
                Unbrs =[neighcol_u(1:(2*b+1)*b+b);neighcol_u((2*b+1)*b+b+2:end)];
                Vnbrs =[neighcol_v(1:(2*b+1)*b+b);neighcol_v((2*b+1)*b+b+2:end)];
                median_Unbrs(j,i)=nanmedian(Unbrs(:));
                median_Vnbrs(j,i)=nanmedian(Vnbrs(:));
                median_u(j,i)=nanmedian(u(:));
                median_v(j,i)=nanmedian(v(:));
                fluct_u(j,i)=delx(j,i)-median_Unbrs(j,i);
                fluct_v(j,i)=dely(j,i)-median_Vnbrs(j,i);
                res_u=Unbrs-median_Unbrs(j,i);
                res_v=Vnbrs-median_Vnbrs(j,i);
                median_res_u=nanmedian(abs(res_u(:)));
                median_res_v=nanmedian(abs(res_v(:)));
                normfluct_u(j,i)=abs(fluct_u(j,i)/(median_res_u+epsilon));
                normfluct_v(j,i)=abs(fluct_v(j,i)/(median_res_v+epsilon));
            end
        end
    end
    reject_median=(sqrt(normfluct_u.^2+normfluct_v.^2)>thresh);
    delx(reject_median==1)=NaN;
    dely(reject_median==1)=NaN;
    delx(reject_median==1)=median_u(reject_median==1);
    dely(reject_median==1)=median_v(reject_median==1);
    delx = delx( b+1:(end-b), b+1:(end-b) );
    dely = dely( b+1:(end-b), b+1:(end-b) );

    %Global fiter - reject displacement more than stdthresh(5) times the standard deviation (of the whole image)
    % [J,I]=size(delx);
    % stdthresh=5;
    % meanu=nanmean(delx(:));
    % meanv=nanmean(dely(:));
    % std2u=nanstd(reshape(delx,J*I,1));
    % std2v=nanstd(reshape(dely,J*I,1));
    % minvalu=meanu-stdthresh*std2u;
    % maxvalu=meanu+stdthresh*std2u;
    % minvalv=meanv-stdthresh*std2v;
    % maxvalv=meanv+stdthresh*std2v;
    % reject_std=(delx<minvalu | delx>maxvalu | dely<minvalv | dely>maxvalv);
    % delx(reject_std==1)=NaN;
    % dely(reject_std==1)=NaN;
    
    %
    % %outlier removal percent estimate
    % reject_total=sum(sum(reject_median==1))+sum(sum(reject_std==1));
    % reject_total_percent=reject_total/(I*J)*100;
    
    %delx=(smoothn(delx,0,'robust'));
    %dely=(smoothn(dely,0,'robust'));
    %dcor=(smoothn(dcor,0,'robust'));
    
    %outlier interpolation
    delx=(smoothn(delx,1/lvl/10,'robust'));
    dely=(smoothn(dely,1/lvl/10,'robust'));
    %dcor=(smoothn(dcor,1/lvl/10,'robust'));
        
    %interpolation for next level
    %grid at next level
    IW2 = IntrWndw(lvl+1); %window size (interogation window)
    GS2 = GrdSpc(lvl+1);
    x2 = IW2/2:GS2:(w-IW2/2);
    y2 = IW2/2:GS2:(h-IW2/2);
    
    %interpolation for next level
    INTdelx = (interp2(x, y, delx, x2, y2', '*spline'));
    INTdely = (interp2(x, y, dely, x2, y2', '*spline'));
    INTdcor = (interp2(x, y, dcor, x2, y2', '*nearest')); %neareast better?
    INTdcor(isnan(INTdcor))=0;
    INTdwhatlevel = (interp2(x, y, dwhatlevel, x2, y2', '*nearest')); %neareast better?
    INTdwhatlevel(isnan(INTdwhatlevel))=0;
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
        VALID_BLOCK=(~isnan(mask(r, c))); %if the center of window is not NAN, it's valid
        gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
        gdely_round = round(INTdely(bxCNTr, bxCNTc));
        
        if VALID_BLOCK %if the PIV can be done (not masked)
            MASK(bxCNTr, bxCNTc)=1;
            if ( (bdryT + gdely_round )>=1 && (bdryB + gdely_round)<=h && ...
                    (bdryL + gdelx_round )>=1 && (bdryR + gdelx_round)<=w ) % check that flow points within the boundary of 2nd image
                
                % sub-window
                bxA = IM1(bdryT:bdryB, bdryL:bdryR);
                bxB = IM2( (bdryT:bdryB) + gdely_round, (bdryL:bdryR) + gdelx_round );
                
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
                Xcorr=fftshift(real(ifft2(conj(fft2(bxAmm)).*fft2(bxBmm))))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); %cross correlation
                [PHpky, PHpkx] = find( PHcorr == max(max(PHcorr)) ); %find max in phase correlation
                [Xpky, Xpkx] = find( Xcorr == max(max(Xcorr)) );%find max in cross correlation
                
                if ( length(PHpky)==1 && length(PHpkx)==1 && length(Xpky)==1 && length(Xpkx)==1 && Xpky==PHpky && Xpkx==PHpkx  ) %if local calculation OK. (both phase and cross correlation agree with single peak)
                    ldelx =  Xpkx - IW/2 - 1; %local velocity calculated here
                    ldely = Xpky - IW/2 - 1;
                    dcor(bxCNTr, bxCNTc) = Xcorr(Xpky, Xpkx); %local correlation calculated here
                    corrmtx(bxCNTr, bxCNTc)= Xcorr(Xpky, Xpkx);
                    if (dcor(bxCNTr, bxCNTc)>0 && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
                        %                 if (dcor(bxCNTr, bxCNTc)>INTdcor(bxCNTr, bxCNTc)-0.5/lvl && abs(ldelx)<IW/2 && abs(ldely)<IW/2 ) %is this correlation better than before
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
                                delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX;
                                dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
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
                    end
                    
                else %no local calculation possible -> keep global velocity
                    %total velocity - global
                    delx(bxCNTr, bxCNTc) = gdelx;
                    dely(bxCNTr, bxCNTc) = gdely;
                    dcor(bxCNTr, bxCNTc) = gcor;
                    dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel;
                end %local calculation possibble
                
                
            else %if current boxes aren't within image boundaries - no velocity measured at level
                %total velocity - global
                delx(bxCNTr, bxCNTc) = gdelx;
                dely(bxCNTr, bxCNTc) = gdely;
                dcor(bxCNTr, bxCNTc) = gcor;
                dwhatlevel(bxCNTr, bxCNTc) = gwhatlevel;
            end %end of boundary check if statement
        end %end VALID BLOCK CHECK
        bxCNTr = bxCNTr + 1;
        
    end    %end of row for loop
    
    bxCNTc = bxCNTc + 1;
    
end  %end of column for loop


%%
% %LOCAL LEVEL FILTERS
% %filter & smooth
%

%Global filter - reject displacement more than stdthresh(5) times the standard deviation of the whole image
% [J,I]=size(delx);
% stdthresh=5;
% meanu=nanmean(delx(:));
% meanv=nanmean(dely(:));
% std2u=nanstd(reshape(delx,J*I,1));
% std2v=nanstd(reshape(dely,J*I,1));
% minvalu=meanu-stdthresh*std2u;
% maxvalu=meanu+stdthresh*std2u;
% minvalv=meanv-stdthresh*std2v;
% maxvalv=meanv+stdthresh*std2v;
% reject_std=(delx<minvalu | delx>maxvalu | dely<minvalv | dely>maxvalv);
% delx(reject_std==1)=NaN;
% dely(reject_std==1)=NaN;
%
% % % Local filter - Reject bad correlations
% reject_correlation=(dcor<0.5);
% delx(reject_correlation==1)=NaN;
% dely(reject_correlation==1)=NaN;

%Local filter - Dynamic mean filter
for fltCNT = 1:2
    b=1; %3x3 filter; b=2 -> 5x5 filter etc...
    delx=padarray(delx,[b b],'replicate','both');
    dely=padarray(dely,[b b],'replicate','both');
    VALID=~isnan(padarray(MASK,[b b],'replicate','both'));
    reject_dynamic_u=delx*0;
    reject_dynamic_v=dely*0;
    c1 = 0.4;
    c2 = 0.4;
    [J,I]=size(delx);
    for i=1+b:I-b
        for j=1+b:J-b
            if VALID(j,i)
                u=delx(j-b:j+b,i-b:i+b);
                v=dely(j-b:j+b,i-b:i+b);
                neighcol_u=u(:);
                neighcol_v=v(:);
                Unbrs =[neighcol_u(1:(2*b+1)*b+b);neighcol_u((2*b+1)*b+b+2:end)];
                Vnbrs =[neighcol_v(1:(2*b+1)*b+b);neighcol_v((2*b+1)*b+b+2:end)];
                if abs( nanmean(Unbrs) - delx(j,i) ) > ( c1 + c2*nanstd(Unbrs) )
                    %if abs( nansum([nanmean(Unbrs) -delx(j,i)] )) > ( c1 + c2*nanstd(Unbrs) )
                    delx(j,i) = nanmedian( neighcol_u );
                    reject_dynamic_u(j,i) = 1;
                end
                if abs( nanmean(Vnbrs) - dely(j,i) ) > ( c1 + c2*nanstd(Vnbrs) )
                    %if abs( nansum([nanmean(Vnbrs) -dely(j,i)] )) > ( c1 + c2*nanstd(Vnbrs) )
                    dely(j,i) = nanmedian( neighcol_v );
                    reject_dynamic_v(j,i) = 1;
                end
            end %VALID
        end    %rows
    end    %columns
    delx = delx( b+1:(end-b), b+1:(end-b) );
    dely = dely( b+1:(end-b), b+1:(end-b) );
    reject_dynamic_u = reject_dynamic_u( b+1:(end-b), b+1:(end-b) );
    reject_dynamic_v = reject_dynamic_v( b+1:(end-b), b+1:(end-b) );
end    %filter count


%Local filter - reject displacement locally where the fluctuation is more than twice the local std
b=1; %3x3 filter; b=2 -> 5x5 filter etc...
delx=padarray(delx,[b b],'replicate','both');
dely=padarray(dely,[b b],'replicate','both');
VALID=~isnan(padarray(MASK,[b b],'replicate','both'));
[J,I]=size(delx);
epsilon=0.1;
%epsilon=0.1*(1/lvl);
thresh=2;
median_Unbrs=NaN(J,I);
median_Vnbrs=NaN(J,I);
median_u=NaN(J,I);
median_v=NaN(J,I);
fluct_u=NaN(J,I);
fluct_v=NaN(J,I);
normfluct_u=NaN(J,I);
normfluct_v=NaN(J,I);
for i=1+b:I-b
    for j=1+b:J-b
        if VALID(j,i)
            u=delx(j-b:j+b,i-b:i+b);
            v=dely(j-b:j+b,i-b:i+b);
            neighcol_u=u(:);
            neighcol_v=v(:);
            Unbrs =[neighcol_u(1:(2*b+1)*b+b);neighcol_u((2*b+1)*b+b+2:end)];
            Vnbrs =[neighcol_v(1:(2*b+1)*b+b);neighcol_v((2*b+1)*b+b+2:end)];
            median_Unbrs(j,i)=nanmedian(Unbrs(:));
            median_Vnbrs(j,i)=nanmedian(Vnbrs(:));
            median_u(j,i)=nanmedian(u(:));
            median_v(j,i)=nanmedian(v(:));
            fluct_u(j,i)=delx(j,i)-median_Unbrs(j,i);
            fluct_v(j,i)=dely(j,i)-median_Vnbrs(j,i);
            res_u=Unbrs-median_Unbrs(j,i);
            res_v=Vnbrs-median_Vnbrs(j,i);
            median_res_u=nanmedian(abs(res_u(:)));
            median_res_v=nanmedian(abs(res_v(:)));
            normfluct_u(j,i)=abs(fluct_u(j,i)/(median_res_u+epsilon));
            normfluct_v(j,i)=abs(fluct_v(j,i)/(median_res_v+epsilon));
        end
    end
end
reject_median=(sqrt(normfluct_u.^2+normfluct_v.^2)>thresh);
delx(reject_median==1)=median_u(reject_median==1);
dely(reject_median==1)=median_v(reject_median==1);
%delx(reject_median==1)=NaN;
%dely(reject_median==1)=NaN;
delx = delx( b+1:(end-b), b+1:(end-b) );
dely = dely( b+1:(end-b), b+1:(end-b) );

%outlier removal percent estimate
%reject_total=sum(sum(reject_median==1))+sum(sum(reject_std==1))+sum(sum(reject_size==1)); %+sum(sum(reject_dynamic==1)) +sum(sum(reject_correlation==1))
%reject_total_percent=reject_total/(I*J)*100;

%outlier interpolation
compVel.delta_x=smoothn(delx,0,'robust'); % displacement in pixels
compVel.delta_z=smoothn(dely,0,'robust'); % displacement in pixels

%even smoother
% delx_ints=smoothn(delx,'robust');
% dely_ints=smoothn(dely,'robust');

compVel.cor_mtx = corrmtx;
compVel.level_mtx = dwhatlevel;
%
compVel.mask = MASK;
%

