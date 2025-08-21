function [CompVel] =  ComputeVelocities_Quick_NoFilt_Deform_Water(PIV1, PIV2, Mask1, Mask2, IntrWndw, GrdSpc)
% The function ComputeVelocities calculates velocities in the air above the
% water (ones in mask), from particle displacement in fused images. This is
% written by Fabrice Veron, and modified later  by Marc Buckley on Novemebr
% 14, 2013.
% 
%    IM1 is fused image a, IM2 is fused image b, Mask is the mask of ones &
%    zeros, IntrWndw is the interrogation  window, and GrdSpc is the number
%    of pixel that interrogation window moves forward. Moreover, delta_x is
%    the horizontal displacements and delta_z is the vertical displacement.
%    The suggested numbers for the interrogation window and grid spaces are
%    as follows and we usually overlap by fifty percent for the grid space:
%    IntrWndw = [ 64 32 16 8];
%    GrdSpc = [ 32 16 8 4];


PIV1 = double(PIV1) .* Mask1;
PIV2 = double(PIV2) .* Mask2;
PIV1(isnan(PIV1)) = nanmean(nanmean(PIV1));
PIV2(isnan(PIV2)) = nanmean(nanmean(PIV2));
IM1_D=PIV1;
IM2_D=PIV2;
% Mask should contain NaNs where there is no possible velocity calculation, 
% ones elsewhere.

[h, w] = size(PIV1); % Image height and width
[X1,Y1] = meshgrid([1:w], [1:h]);

number_of_levels = length(IntrWndw); % The first four, five, six levels are 
% called the "Global Level" and the last level is called the "Local Level".

%%
% ----------------------- GOBAL LEVEL CALCULATIONS ---------------------- %
for lvl = 1:number_of_levels-1 % First level
    
    IW = IntrWndw(lvl); % Interrogation window size
    GS = GrdSpc(lvl); % Grid spacing size
    
    x = IW/2:GS:(w-IW/2); % x-coordinate of center of interrogation windows
    y = IW/2:GS:(h-IW/2); % y-coordinate of center of interrogation windows
    
    bxsNh = floor(1 + (h - IW)/GS); % No. of interrogation blocks in height
    bxsNw = floor(1 + (w - IW)/GS); % No. of interrogation blocks in width
    
    delx = NaN(bxsNh, bxsNw); % Displacement in x - total velocity in x
    dely = NaN(bxsNh, bxsNw); % Displacement in y - total velocity in y
    dcor = NaN(bxsNh, bxsNw); % Correlation
    % i.e. where there is data .
    
    
    
    if lvl == 1 % Initialize the global motion and correlation if it is the
    % first level , otherwise it is passed from previous  level after being
    % interpolated to proper grid.
        pdelx = zeros(bxsNh, bxsNw);
        pdely = zeros(bxsNh, bxsNw);
    end
    

    bxCNTc = 1; % Column counter initialization
        
    for c = x % Loop in column: x-coordinate of interrogation window
        
        bxCNTr = 1; % Row counter initialization
        
        bdryL = c - IW/2 + 1; % Left and right boundaries of sub-window
        bdryR = c + IW/2;
        
        for r = y % Loop in row: y-coordinate of interrogation window
            
            bdryT = r - IW/2 + 1; % Top and bottom boundaries of sub-window
            bdryB = r + IW/2;
            
            % Global velocity and correlation in this block , which is from
            % previous level, from what we have in memory.
            
            VALID_BLOCK = (~isnan(Mask1(r,c))); % If center of window is NaN
            
            % velocity to displace the interrogation box in the 2nd image
            
            % How far forward or backward the next window is gona be pushed
                 
            
            if VALID_BLOCK % Whether the PIV can be done (not masked)
                try
                    
                    % Sub-window blocks
                    bxA = IM1_D( (bdryT:bdryB), (bdryL:bdryR));
                    bxB = IM2_D( (bdryT:bdryB), (bdryL:bdryR));
                    
                    % Sub-window demeaned, windowed, and FFTs
                    bxAmm = bxA - mean(bxA(:));
                    bxBmm = bxB - mean(bxB(:));
                    
                    fftA = fft2(bxAmm);
                    fftB = fft2(bxBmm);
                    fftCorr = fftB .* conj(fftA);
                    
                                         
                    Xcorr = fftshift(real(ifft2(fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); % Cross Correlation
                    [Xpky,Xpkx] = find( Xcorr == max(max(Xcorr)) ); % Find max in the cross correlation
                    
                    ldelx =  Xpkx - IW/2 - 1; % Local velocity calculated here
                    ldely =  Xpky - IW/2 - 1;
                    
                    % 3 Point Gaussian Interpolation
                     T = log(Xcorr (Xpky-1:Xpky+1, Xpkx-1:Xpkx+1) );
                     t = T(:,2);
                     SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                     t = T(2,:);
                     SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));

           

                     if (isreal([SubpixelY SubpixelX]) && (SubpixelY<1) && (SubpixelX<1) &&  ldelx< IW/2 &&  ldely< IW/2)
                     %if ( length(Xpky)==1 && length(Xpkx)==1 &&  ldelx< IW/2 &&  ldely< IW/2)
                        % Velocity is round of  previous guess (do not keep 
                        % previous subpixel estimate) + local + subpix
                        delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc) + ldelx + SubpixelX;                           
                        dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc) + ldely + SubpixelY;
                        %GLOBAL(bxCNTr, bxCNTc) = 1;
                        dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                   
                    end
                    
                end
            
            end % End of valid block check
            
            bxCNTr = bxCNTr + 1;
            
        end % End of row for loop
        
        bxCNTc = bxCNTc + 1;
        
    end % End of column for loop
    
    % ---------------------- GLOBAL LEVEL FILTERS ----------------------- %
    % Filter and smooth, reject displacement more than 1/2 the window size
    
    % Outlier interpolation
    delx = (smoothn(delx,'robust'));
    dely = (smoothn(dely,'robust'));
    
    
    [X,Y] = meshgrid(x, y);
    %Interpolation of velocity to image resolution
    U1 = interp2(X,Y,delx,X1,Y1,'*spline'); 
    V1 = interp2(X,Y,dely,X1,Y1,'*spline');
    
    %Warping both images according to velocity (centered difference)
    IM1_D= interp2(1:size(PIV1,2),(1:size(PIV1,1))',PIV1,X1,Y1,'*linear');
    IM2_D= interp2(1:size(PIV2,2),(1:size(PIV2,1))',PIV2,X1+U1,Y1+V1,'*linear');
    
    %Warping second image according to velocity
    %IM2_D= interp2(1:size(IM2,2),(1:size(IM2,1))',IM2,X1+U1,Y1+V1,'*linear');
    %IM1_D=IM1;
    
    % Interpolation velocity field for the next level
    % Grid at next level
    IW2 = IntrWndw(lvl + 1); % Window size (interogation window)
    GS2 = GrdSpc(lvl + 1);
    x2 = IW2/2:GS2:(w - IW2/2);
    y2 = IW2/2:GS2:(h - IW2/2);
    % Interpolation was already done (downsampling from image resolution)
    pdelx = U1(y2,x2);
    pdely = V1(y2,x2);

end

% -------------------------- END GLOBAL MOTION -------------------------- %

%%
% -------------------- LOCAL (LAST) LEVEL CALCULATION ------------------- %

lvl = number_of_levels;

IW = IntrWndw(lvl); % Interrogation window size
GS = GrdSpc(lvl); % Grid spacing size

x = IW/2:GS:(w - IW/2); % x-coordinate of center of interrogation windows
y = IW/2:GS:(h - IW/2); % x-coordinate of center of interrogation windows

bxsNh = floor(1 + (h - IW)/GS); % Number of interrogation blocks in height
bxsNw = floor(1 + (w - IW)/GS); % Number of interrogation blocks in width

delx = pdelx;% NaN(bxsNh, bxsNw); % Total velocity in x - accumulates from previous levels with INTdel and gdel
dely = pdely;%NaN(bxsNh, bxsNw); % Total velocity in y - accumulates from previous levels with INTdel and gdel
dcor = NaN(bxsNh, bxsNw); % Correlation - accumulates from previous levels with INTcor and gcor
MASK = NaN(bxsNh, bxsNw); % Shows where the PIV calculation is valid , i.e. where there's data

bxCNTc = 1; % Column counter initialization

for c = x % Loop in column: x-coordinate of center of interrogation window
    
    bxCNTr = 1; % Row counter initialization
    
    bdryL = c - IW/2 + 1;  % Left and right boundaries of sub-window
    bdryR = c + IW/2;
    
    for r = y % Loop in row: y-coordinate of center of interrogation window
        
        bdryT = r - IW/2 + 1; % Top and bottom boundaries of sub-window
        bdryB = r + IW/2;
        
        % Global velocity & correlation in this block (from previous level)
        VALID_BLOCK = (~isnan(Mask1(r,c))); % If no NANs at all in subwindow
       
        % velocity to displace the interrogation box in the 2nd image
        
       
        
        if VALID_BLOCK % Whether the PIV can be done (not masked)
            MASK(bxCNTr, bxCNTc) = 1;
            try
                
                % Sub-window blocks
                bxA = IM1_D( (bdryT:bdryB), (bdryL:bdryR));
                bxB = IM2_D( (bdryT:bdryB), (bdryL:bdryR));
                
                % Sub-window demeaned, windowed, and FFTs
                bxAmm = bxA - mean(bxA(:));
                bxBmm = bxB - mean(bxB(:));
                
                fftA = fft2(bxAmm);
                fftB = fft2(bxBmm);
                fftCorr = fftB .* conj(fftA);
                
                % PHcorr = fftshift(abs(ifft2(fftCorr./(abs(fftCorr)+eps)))); % Phase correlation
                % [PHpky,PHpkx] = find(PHcorr==max(max(PHcorr))); % Find max in phase correlation
                
                Xcorr = fftshift(real(ifft2(fftCorr)))./sqrt(sum(sum(bxAmm.^2)))./sqrt(sum(sum(bxBmm.^2))); % Cross correlation
                [Xpky,Xpkx] = find( Xcorr == max(max(Xcorr)) ); % Find max in cross correlation
                
                ldelx = Xpkx - IW/2 - 1; % Local velocity calculated here
                ldely = Xpky - IW/2 - 1;

                % 3 Point Gaussian Interpolation
                 T = log(Xcorr (Xpky-1:Xpky+1, Xpkx-1:Xpkx+1)); 
                 t = T(:,2);
                 SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                 t = T(2,:);
                 SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));

                
                if (isreal([SubpixelY SubpixelX]) && (SubpixelY<1) && (SubpixelX<1) )
                    % Velocity is round of previous guess (do not keep pre-
                    % vious subpixel estimate) + local + subpix
                    delx(bxCNTr, bxCNTc) = pdelx(bxCNTr, bxCNTc) + ldelx + SubpixelX;
                    dely(bxCNTr, bxCNTc) = pdely(bxCNTr, bxCNTc) + ldely + SubpixelY;
                    dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                
                end
          
            end
            
        end % End of valid block check
        
        bxCNTr = bxCNTr + 1;
        
    end % End of row for loop
    
    bxCNTc = bxCNTc + 1;
    
end % End of column for loop

[X,Y] = meshgrid(x, y);

U1 = interp2(X,Y,delx,X1,Y1,'*spline'); 
V1 = interp2(X,Y,dely,X1,Y1,'*spline');

% --------------------- END LOCAL LEVEL CALCULATION --------------------- %
%%

CompVel.INTdelx = pdelx;
CompVel.INTdelz = - pdely;
CompVel.dcor = dcor;
CompVel.Mask = MASK;

CompVel.delta_x = delx; % Final Result: vectors going downwind (from left to right) are positive
CompVel.delta_z = - dely; % Final Result: vectors going away from surface (upward) are positive

CompVel.delta_x1 = U1; %Final Result at resolution of the original image
CompVel.delta_z1 = -V1; %Final Result at resolution of the original image

CompVel.xPIV = x; % 0 is upper left corner of image (upwind, away from surface)
CompVel.zPIV = y; % 0 is upper left corner of image (upwind, away from surface)
CompVel.GS = GrdSpc(end); % Final grid spacing
CompVel.IW = IntrWndw(end);

end

