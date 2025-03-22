function [CompVel] =  ComputeVelocities_Quick_NoFilt_Water_InitVelOrb(IM1, IM2, Mask1, Mask2, Mask, IntrWndw, GrdSpc,Uorb_W1,Vorb_W1)
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


IM1 = double(IM1) .* Mask1;
IM2 = double(IM2) .* Mask2;
IM1(isnan(IM1)) = nanmean(nanmean(IM1));
IM2(isnan(IM2)) = nanmean(nanmean(IM2));
% Mask should contain NaNs where there is no possible velocity calculation,
% ones elsewhere.

[h, w] = size(IM1); % Image height and width

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
    GLOBAL = NaN(bxsNh, bxsNw); % Shows where the GLOBAL Motion is calculated at this level only
    MASK = NaN(bxsNh, bxsNw); % Shows where the PIV calculation is valid,
    % i.e. where there is data .



    if lvl == 1 % Initialize the global motion and correlation if it is the
        % first level , otherwise it is passed from previous  level after being
        % interpolated to proper grid.
        for i = 1:length(x)
            for ii = 1:length(y)
                Uguess(ii,i) = mean(mean(Uorb_W1(y(ii)-IW/2+1:y(ii)+IW/2-1,x(i)-IW/2+1:x(i)+IW/2-1),'omitnan'),'omitnan');
                Vguess(ii,i) = mean(mean(Vorb_W1(y(ii)-IW/2+1:y(ii)+IW/2-1,x(i)-IW/2+1:x(i)+IW/2-1),'omitnan'),'omitnan');
            end
        end
        INTdelx = Uguess; %zeros(bxsNh, bxsNw);
        INTdely = Vguess; %zeros(bxsNh, bxsNw);
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

            VALID_BLOCK = (~isnan(Mask(r,c))); % If center of window is NaN
            gdelx_round = round(INTdelx (bxCNTr,bxCNTc) );
            gdely_round = round(INTdely (bxCNTr,bxCNTc) ); % Roundup global
            % velocity to displace the interrogation box in the 2nd image

            % How far forward or backward the next window is gona be pushed
            forwardx = ceil(gdelx_round);
            forwardy = ceil(gdely_round);



            if VALID_BLOCK % Whether the PIV can be done (not masked)
                MASK(bxCNTr, bxCNTc) = 1;
                try

                    % Sub-window blocks
                    bxA = IM1( (bdryT:bdryB) , (bdryL:bdryR) );
                    bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );

                    % Sub-window demeaned, windowed, and FFTs
                    bxAmm = bxA - mean(bxA(:));
                    bxBmm = bxB - mean(bxB(:));

                    fftA = fft2(bxAmm);
                    fftB = fft2(bxBmm);
                    fftCorr = fftB .* conj(fftA);

                    %PHcorr = fftshift(abs(ifft2(fftCorr./(abs(fftCorr)+eps)))); % Phase Correlation
                    %[PHpky,PHpkx] = find(PHcorr==max(max(PHcorr))); % Find max in phase correlation
                    %[Xpky,Xpkx] = find( PHcorr == max(max(PHcorr)) );

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



                    if (isreal(SubpixelY) && isreal(SubpixelX) && (SubpixelY<1) && (SubpixelX<1) &&  ldelx< IW/2 &&  ldely< IW/2)
                        %if ( length(Xpky)==1 && length(Xpkx)==1 &&  ldelx< IW/2 &&  ldely< IW/2)
                        % Velocity is round of  previous guess (do not keep
                        % previous subpixel estimate) + local + subpix
                        delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX;
                        dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                        %GLOBAL(bxCNTr, bxCNTc) = 1;
                        %dcor(bxCNTr, bxCNTc) = max(max(Xcorr));

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


    % Interpolation for the next level
    % Grid at next level
    IW2 = IntrWndw(lvl + 1); % Window size (interogation window)
    GS2 = GrdSpc(lvl + 1);
    x2 = IW2/2:GS2:(w - IW2/2);
    y2 = IW2/2:GS2:(h - IW2/2);

    % Interpolation for the next level
    INTdelx = (interp2(x, y', delx, x2, y2', '*spline'));
    INTdely = (interp2(x, y', dely, x2, y2', '*spline'));

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

delx = NaN(bxsNh, bxsNw); % Total velocity in x - accumulates from previous levels with INTdel and gdel
dely = NaN(bxsNh, bxsNw); % Total velocity in y - accumulates from previous levels with INTdel and gdel
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
        VALID_BLOCK = (~isnan(Mask(r,c))); % If no NANs at all in subwindow
        gdelx_round = round(INTdelx(bxCNTr, bxCNTc));
        gdely_round = round(INTdely(bxCNTr, bxCNTc)); % Round-up the global
        % velocity to displace the interrogation box in the 2nd image

        forwardx = ceil(gdelx_round);
        forwardy = ceil(gdely_round);

        ldelx = nan; % Local velocity calculated here
        ldely = nan;

        if VALID_BLOCK % Whether the PIV can be done (not masked)
            MASK(bxCNTr, bxCNTc) = 1;
            try

                % Sub-window blocks
                bxA = IM1( (bdryT:bdryB) , (bdryL:bdryR) );
                bxB = IM2( (bdryT:bdryB) + forwardy, (bdryL:bdryR) + forwardx );

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
                
                if length(Xpkx) == 1
                    ldelx = Xpkx - IW/2 - 1; % Local velocity calculated here
                    ldely = Xpky - IW/2 - 1;
                end

                % T = (Xcorr (Xpky-1:Xpky+1, Xpkx-1:Xpkx+1));
                % [X,Y] = meshgrid(-1:1,-1:1);
                % [Xq,Yq] = meshgrid(-1:0.01:1,-1:0.01:1);
                % Vq = interp2(X,Y,T,Xq,Yq,'*spline');
                % xx=[-1:0.01:1]; SubpixelY=xx(i); SubpixelX=xx(j);

                % 3 Point Gaussian Interpolation
                T = log(Xcorr (Xpky-1:Xpky+1, Xpkx-1:Xpkx+1));

                t = T(:,2);
                SubpixelY = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                t = T(2,:);
                SubpixelX = (1/2)*(t(3) - t(1))/(2*t(2) - t(1) - t(3));
                %         for i = -1:1
                %             for j = -1:1
                %
                %                 c10(j+2,i+2)=i*T (i+2,j+2);
                %                 c01(j+2,i+2) = j*T (i+2,j+2);
                %                 c11(j+2,i+2) = i*j*T (i+2,j+2);
                %                 c20(j+2,i+2) = (3*i^2-2)*T (i+2,j+2);
                %                 c02(j+2,i+2) = (3*j^2-2)*T (i+2,j+2);
                %             end
                %         end
                %         c10 = (1/6)*sum(sum(c10));
                %         c01 = (1/6)*sum(sum(c01));
                %         c11 = (1/4)*sum(sum(c11));
                %         c20 = (1/6)*sum(sum(c20));
                %         c02 = (1/6)*sum(sum(c02));
                %         SubpixelX = squeeze((c11.*c01-2*c10.*c02)./(4*c20.*c02-c11.^2));
                %         SubpixelY = squeeze((c11.*c10-2*c01.*c20)./(4*c20.*c02-c11.^2));
                
            end %previously, end of try was after the if statement below
            
            if VALID_BLOCK
                if (isreal(SubpixelY) && isreal(SubpixelX) && (SubpixelY<1) && (SubpixelX<1) &&  ldelx< IW/2 &&  ldely< IW/2)
                    % Velocity is round of previous guess (do not keep pre-
                    % vious subpixel estimate) + local + subpix
                    delx(bxCNTr, bxCNTc) = gdelx_round + ldelx + SubpixelX;
                    dely(bxCNTr, bxCNTc) = gdely_round + ldely + SubpixelY;
                    dcor(bxCNTr, bxCNTc) = max(max(Xcorr));

                    %Elseif and else added by Andy to eliminate NaNs in
                    %final output.
                elseif ~isnan(ldelx)&& ldelx < IW/2 &&  ldely < IW/2
                    delx(bxCNTr, bxCNTc) = gdelx_round + ldelx;
                    dely(bxCNTr, bxCNTc) = gdely_round + ldely;
                    dcor(bxCNTr, bxCNTc) = max(max(Xcorr));
                else
                    delx(bxCNTr, bxCNTc) = gdelx_round;
                    dely(bxCNTr, bxCNTc) = gdely_round;
                end
            end

            

        end % End of valid block check

        bxCNTr = bxCNTr + 1;

    end % End of row for loop

    bxCNTc = bxCNTc + 1;

end % End of column for loop

% --------------------- END LOCAL LEVEL CALCULATION --------------------- %
%%

CompVel.INTdelx = INTdelx;
CompVel.INTdelz = - INTdely;
CompVel.dcor = dcor;
CompVel.Mask = MASK;

CompVel.delta_x = delx; % Final Result: vectors going downwind (from left to right) are positive
CompVel.delta_z = - dely; % Final Result: vectors going away from surface (upward) are positive

CompVel.xPIV = x; % 0 is upper left corner of image (upwind, away from surface)
CompVel.zPIV = y; % 0 is upper left corner of image (upwind, away from surface)
CompVel.GS = GrdSpc(end); % Final grid spacing

end

