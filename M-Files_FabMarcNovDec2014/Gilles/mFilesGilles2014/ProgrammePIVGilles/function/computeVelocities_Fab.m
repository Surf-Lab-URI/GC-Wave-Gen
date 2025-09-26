
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
% 11/14/2013

%%
% IM1=IM_a.*mask; 
% IM2=IM_b.*mask;
%IntrWndw=[128 64 32 16 8];
% GrdSpc=[64 32 16 8 4];
%IntrWndw=[32 16 8];
% GrdSpc=[16 8 4];

function compVel =  computeVelocities_Fab(IM1, IM2, mask1, mask2, IntrWndw, GrdSpc)
% [delx_int dely_int, MASK,cor_mtx,level_mtx]
IM1=double(IM1).*mask1;
IM2=double(IM2).*mask2;
IM1(isnan(IM1))=nanmean(nanmean(IM1));
IM2(isnan(IM2))=nanmean(nanmean(IM2));
% mask should contain NaN where there's no possible velocity calculation, 1s elsewhere

[h, w] = size(IM1); %image height and width
% IntrWndw=[32 16 8];
% GrdSpc=[16 8 4]; 
% IntrWndw=[128 64 32 16 8];
% GrdSpc=[64 32 16 8 4];
number_of_levels= length(IntrWndw);
corr_threshold=0.4;
corr_min=0.25;
mask=mask1; %forward diff

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
    IM_MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
    
    
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
            
           
            IM_MASK(bxCNTr, bxCNTc)= (~isnan(mask1(r, c))); %if the center of window is not NAN, it's valid
            %VALID_BLOCK = ~isnan(sum(sum(mask(bdryT:bdryB,bdryL:bdryR)))); %if no NANs at all in subwindow
            VALID_BLOCK = (~isnan(mask(r, c))); %if no NANs at all in subwindow
            gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
            gdely_round = round(INTdely(bxCNTr, bxCNTc));
            
            forwardx=ceil(gdelx_round/2)*2; %forward diff
            forwardy=ceil(gdely_round/2)*2;%forward diff 
            backx=floor(gdelx_round/2)*0; %forward diff
            backy=floor(gdely_round/2)*0; %forward diff
            
            if VALID_BLOCK %if the PIV can be done (not masked)
                MASK(bxCNTr, bxCNTc)=1;
                try
                    % sub-window blocks
                    
                    bxA = IM1( (bdryT:bdryB) - backy, (bdryL:bdryR) - backx );
                    bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );
                    %hanning sub-window
                    N=size(bxA);
                    Nx=N(2);
                    Ny=N(1);
                    Hanning2d=0.5*(1-cos(2*pi*(0:Ny-1)'/(Ny-1)))*0.5*(1-cos(2*pi*(0:Nx-1)/(Nx-1)));
                    
                    % sub-window demeaned, windowed, and FFTs
                    bxAmm = bxA-mean(bxA(:));
                    bxBmm = bxB-mean(bxB(:));
                     %bxAmm = bxAmm.*Hanning2d;
                     %bxBmm = bxBmm.*Hanning2d;
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
                                if (isreal(SubpixelY+SubpixelX))
                                
                                delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                                GLOBAL(bxCNTr, bxCNTc) = 1;
                                dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                                else
                                delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                dely(bxCNTr, bxCNTc) = NaN;
                                end
                                
                %catch
                 %               delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                 %               dely(bxCNTr, bxCNTc) = NaN;
                end
                    
                  
                    
                    
                
            end %end VALID BLOCK CHECK
            bxCNTr = bxCNTr + 1;
            
        end    %end of row for loop
        
        bxCNTc = bxCNTc + 1;
        
    end  %end of column for loop
    
    % ---------- GLOBAL LEVEL FILTERS -------------------
    %filter & smooth
    %reject displacement more than 1/2 the window size
%    corr_threshold_at_level=corr_threshold-IntrWndw(lvl)/IntrWndw(1)*(corr_threshold-corr_min);
    %i=isnan(MASK);
    %delx(i)=0; dely(i)=0;
    
%   delx(dcor<corr_threshold_at_level)=NaN;
%   dely(dcor<corr_threshold_at_level)=NaN;
   
%    b=1; %3x3 filter; b=2 -> 5x5 filter etc...
%     delx=padarray(delx,[b b],'replicate','both');
%     dely=padarray(dely,[b b],'replicate','both');
%     valid=~isnan(padarray(MASK,[b b],'replicate','both'));
%     [J,I]=size(delx);
%     epsilon=0.1*(1/lvl); %becomes more and more stringent with level.
%     epsilon=0.1; 
%     thresh=5;
%     median_Unbrs=NaN(J,I);
%     median_Vnbrs=NaN(J,I);
%     median_u=NaN(J,I);
%     median_v=NaN(J,I);
%     fluct_u=NaN(J,I);
%     fluct_v=NaN(J,I);
%     normfluct_u=NaN(J,I);
%     normfluct_v=NaN(J,I);
%     for i=1+b:I-b
%         for j=1+b:J-b
%             if valid(j,i)
%                 u=delx(j-b:j+b,i-b:i+b);
%                 v=dely(j-b:j+b,i-b:i+b);
%                 neighcol_u=u(:);
%                 neighcol_v=v(:);
%                 Unbrs =[neighcol_u(1:(2*b+1)*b+b);neighcol_u((2*b+1)*b+b+2:end)];
%                 Vnbrs =[neighcol_v(1:(2*b+1)*b+b);neighcol_v((2*b+1)*b+b+2:end)];
%                 median_Unbrs(j,i)=nanmedian(Unbrs(:));
%                 median_Vnbrs(j,i)=nanmedian(Vnbrs(:));
%                 median_u(j,i)=nanmedian(u(:));
%                 median_v(j,i)=nanmedian(v(:));
%                 fluct_u(j,i)=delx(j,i)-median_Unbrs(j,i);
%                 fluct_v(j,i)=dely(j,i)-median_Vnbrs(j,i);
%                 res_u=Unbrs-median_Unbrs(j,i);
%                 res_v=Vnbrs-median_Vnbrs(j,i);
%                 median_res_u=nanmedian(abs(res_u(:)));
%                 median_res_v=nanmedian(abs(res_v(:)));
%                 normfluct_u(j,i)=abs(fluct_u(j,i)/(median_res_u+epsilon));
%                 normfluct_v(j,i)=abs(fluct_v(j,i)/(median_res_v+epsilon));
%             end
%         end
%     end
%     reject_median=(sqrt(normfluct_u.^2+normfluct_v.^2)>thresh);
%     delx(reject_median==1)=NaN;
%     dely(reject_median==1)=NaN;
%     delx = delx( b+1:(end-b), b+1:(end-b) );
%     dely = dely( b+1:(end-b), b+1:(end-b) );
%     
    
    
   %delx=inpaint_nans(delx,4);
   %dely=inpaint_nans(dely,4);
corr_threshold_at_level=corr_threshold-IntrWndw(lvl)/IntrWndw(1)*(corr_threshold-corr_min);
delx(dcor<corr_threshold_at_level)=NaN;
dely(dcor<corr_threshold_at_level)=NaN;

    %outlier interpolation
    warning off
    delx=(smoothn(delx,'robust'));
    dely=(smoothn(dely,'robust'));
    
    %outlier interpolation
    %delx=(smoothn(delx,'robust'));
    %dely=(smoothn(dely,'robust'));
    %dcor=(smoothn(dcor,1/lvl/10,'robust'));
    
    %interpolation for next level
    %grid at next level
    IW2 = IntrWndw(lvl+1); %window size (interogation window)
    GS2 = GrdSpc(lvl+1);
    x2 = IW2/2:GS2:(w-IW2/2);
    y2 = IW2/2:GS2:(h-IW2/2);
    
    %interpolation for next level
    INTdelx = (interp2(x, y', delx, x2, y2', '*spline'));
    INTdely = (interp2(x, y', dely, x2, y2', '*spline'));
    
end

delx_global=delx;
dely_global=dely;
mask_global=MASK;

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
MASK = NaN(bxsNh, bxsNw); % shows where the PIV calculation is valid (where there's data).
IM_MASK = NaN(bxsNh, bxsNw); 


bxCNTc = 1; %counter initialization

for c = x %loop in column (x, coordinate in image)
        
        bxCNTr = 1; %counter initialization (row)
        
        bdryL = c - IW/2 + 1;  % left and right boundaries of sub-window
        bdryR = c + IW/2;
        
        for r = y %loop in row (y, coordinate in image)
            
            bdryT = r - IW/2 + 1; % top and bottom boundaries of sub-window
            bdryB = r + IW/2;
            
            %global velocity (and correlation) in this block (from previous level)
            %gdelx = INTdelx(bxCNTr, bxCNTc);
            %gdely = INTdely(bxCNTr, bxCNTc);      
            IM_MASK(bxCNTr, bxCNTc)= (~isnan(mask(r, c))); %if the center of window is not NAN, it's valid
            %VALID_BLOCK = ~isnan(sum(sum(mask(bdryT:bdryB,bdryL:bdryR)))); %if no NANs at all in subwindow
            VALID_BLOCK = (~isnan(mask(r, c))); %if no NANs at all in subwindow
            gdelx_round = round(INTdelx(bxCNTr, bxCNTc)); %round-up global velocity to displace the interrogation box in the 2nd image
            gdely_round = round(INTdely(bxCNTr, bxCNTc));
            
            forwardx=ceil(gdelx_round/2)*2;
            forwardy=ceil(gdely_round/2)*2;
            backx=floor(gdelx_round/2)*0;
            backy=floor(gdely_round/2)*0;
            
            if VALID_BLOCK %if the PIV can be done (not masked)
                MASK(bxCNTr, bxCNTc)=1;
                try
                    % sub-window blocks
                    
                    bxA = IM1( (bdryT:bdryB) - backy, (bdryL:bdryR) - backx );
                    bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );
                    %hanning sub-window
                    N=size(bxA);
                    Nx=N(2);
                    Ny=N(1);
                    Hanning2d=0.5*(1-cos(2*pi*(0:Ny-1)'/(Ny-1)))*0.5*(1-cos(2*pi*(0:Nx-1)/(Nx-1)));
                    
                    % sub-window demeaned, windowed, and FFTs
                    bxAmm = bxA-mean(bxA(:));
                    bxBmm = bxB-mean(bxB(:));
                    %bxAmm = bxAmm.*Hanning2d;                    
                    %bxBmm = bxBmm.*Hanning2d;
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
                                if (isreal(SubpixelY+SubpixelX))
                                delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                                dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                                else
                                delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                dely(bxCNTr, bxCNTc) = NaN;
                                %delx(bxCNTr, bxCNTc) = gdelx_round + ldelx; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                %dely(bxCNTr, bxCNTc) = gdely_round + ldely;
                                end
                  %catch              
                
                                %delx(bxCNTr, bxCNTc) = NaN; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                %dely(bxCNTr, bxCNTc) = NaN;
                                %delx(bxCNTr, bxCNTc) = gdelx; %velocity is round of previous guess (do not keep previous subpixel estimate) +local+subpix
                                %dely(bxCNTr, bxCNTc) = gdely;
                end
                    
                  
                    
                    
                
            end %end VALID BLOCK CHECK
            bxCNTr = bxCNTr + 1;
            
        end    %end of row for loop
        
        bxCNTc = bxCNTc + 1;
        
    end  %end of column for loop
   


%%
% %LOCAL LEVEL FILTERS
% %filter & smooth
%
corr_threshold_at_level=corr_threshold-IntrWndw(lvl)/IntrWndw(1)*(corr_threshold-corr_min);
delx(dcor<corr_threshold_at_level)=NaN;
dely(dcor<corr_threshold_at_level)=NaN;
dcor(dcor<corr_threshold_at_level)=NaN;


rejx=abs(delx-INTdelx);
rejy=abs(dely-INTdely);
delx(rejx>2)=NaN;
dely(rejy>2)=NaN;

%delx(isnan(delx))=INTdelx(isnan(delx));
%dely(isnan(dely))=INTdely(isnan(dely));
    %i=isnan(MASK);
    %delx(i)=0; dely(i)=0;
    %delx=inpaint_nans(delx,1);
    %dely=inpaint_nans(dely,1);

    %Mirror
% for i=1:length(delx)
% m=MASK(:,i);
% l=find(isnan(m),1,'first');
%     for k=1:number_of_levels
%     delx(l+k-1,i)=-delx(l-k,i);
%     dely(l+k-1,i)=-dely(l-k,i);
%     end
% end
% 
% b=1; %3x3 filter; b=2 -> 5x5 filter etc...
%     delx=padarray(delx,[b b],'replicate','both');
%     dely=padarray(dely,[b b],'replicate','both');
%     valid=~isnan(padarray(MASK,[b b],'replicate','both'));
%     [J,I]=size(delx);
%     epsilon=0.1*(1/lvl); %becomes more and more stringent with level.
%     epsilon=0.03;
%     thresh=2;
%     normfluct_u=NaN(J,I);
%     normfluct_v=NaN(J,I);
%     for i=1+b:I-b
%         for j=1+b:J-b
%             if valid(j,i)
%                 u=delx(j-b:j+b,i-b:i+b);
%                 v=dely(j-b:j+b,i-b:i+b);
%                 neighcol_u=u(:);
%                 neighcol_v=v(:);
%                 Unbrs =[neighcol_u(1:(2*b+1)*b+b);neighcol_u((2*b+1)*b+b+2:end)];
%                 Vnbrs =[neighcol_v(1:(2*b+1)*b+b);neighcol_v((2*b+1)*b+b+2:end)];
%                 median_Unbrs=nanmedian(Unbrs(:));
%                 median_Vnbrs=nanmedian(Vnbrs(:));
%                 fluct_u=delx(j,i)-median_Unbrs;
%                 fluct_v=dely(j,i)-median_Vnbrs;
%                 res_u=Unbrs-median_Unbrs;
%                 res_v=Vnbrs-median_Vnbrs;
%                 median_res_u=nanmedian(abs(res_u(:)));
%                 median_res_v=nanmedian(abs(res_v(:)));
%                 normfluct_u(j,i)=abs(fluct_u/(median_res_u+epsilon));
%                 normfluct_v(j,i)=abs(fluct_v/(median_res_v+epsilon));
%             end
%         end
%     end
%     reject_median=(sqrt(normfluct_u.^2+normfluct_v.^2)>thresh);
%     delx(reject_median==1)=NaN;
%     dely(reject_median==1)=NaN;
%     delx = delx( b+1:(end-b), b+1:(end-b) );
%     dely = dely( b+1:(end-b), b+1:(end-b) );
    
    %outlier interpolation
    delx=(smoothn(delx,'robust'));
    dely=(smoothn(dely,'robust'));
 
% Linear extrap    
%  ddx = nan(size(delx));
%  ddy = nan(size(dely));
%  y = 1:size(delx,1);
%  for i=1:size(delx,2)
%      tempx=(delx(:,i));
%      tempy=(dely(:,i));
%      nx=~isnan(tempx);
%      ny=~isnan(tempy);
%      yx=y(nx);
%      yy=y(ny);
%      tempx=tempx(nx);
%      tempy=tempy(ny);
%      try
%          ax_i=interp1(yx,tempx,y,'linear','extrap');
%          ay_i=interp1(yy,tempy,y,'linear','extrap');
%          ddx(:,i)=ax_i;
%          ddy(:,i)=ay_i;
%      catch
%          ddx(:,i)=NaN;
%          ddy(:,i)=NaN;
%      end
%  end
%  
%  delx = ddx; 
%  dely = ddy; 
 

compVel.dcor = dcor;
compVel.mask = MASK;
IM_MASK(IM_MASK==0)=NaN;
compVel.immask = IM_MASK;


compVel.delta_x=delx;
compVel.delta_z=dely;
