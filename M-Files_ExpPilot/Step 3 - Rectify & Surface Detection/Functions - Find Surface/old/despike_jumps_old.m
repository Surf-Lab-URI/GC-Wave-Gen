function surf2_corr = despike_jumps_old(surface2,Jump_thr,MarkInt,Type)

%%%%%% Function used to despike Shoaling waves from sudden jumps in the
%%%%%% surface detection

% Step 0: re-interpolate previously interpolated points
surf2_corr = surface2;
for i = 1:length(MarkInt)
    if length(MarkInt)>1
        if MarkInt(i)>1500 && MarkInt(i)<9000
            VInt = polyval(polyfit(MarkInt(i)-250:MarkInt(i)+250,surf2_corr(MarkInt(i)-250:MarkInt(i)+250),30),MarkInt(i)-250:MarkInt(i)+250);
            surf2_corr(MarkInt(i)-25:MarkInt(i)+25) = VInt(226:end-225);
        end
    end
end
% Step 1: small interval of spikes
[~,locs1] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);
DiffSurf = diff(surf2_corr);
IntLocs = [];
c = 0;
for ii = 1:length(locs1)-1
    if ii == c
        continue
    end
    if DiffSurf(locs1(ii))*DiffSurf(locs1(ii+1))<0 && abs(locs1(ii)-locs1(ii+1))<100
        IntLocs = [IntLocs,locs1(ii):locs1(ii+1)-1];
        c = ii+1;
    end
end
IntLocs = IntLocs+1;
surf2_corr = surface2;
surf2_corr(IntLocs) = NaN;
surf2_corr(601:end) = interp1(601:length(surf2_corr),surf2_corr(601:end),601:length(surf2_corr),'spline','extrap');

% Step 2: bigger interval of spikes
[~,locs2] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);
DiffSurf = diff(surf2_corr);
IntLocs = [];
c = 0;
for ii = 1:length(locs2)-1
    if ii == c
        continue
    end
    if DiffSurf(locs2(ii))*DiffSurf(locs2(ii+1))<0 && abs(locs2(ii)-locs2(ii+1))<200
        IntLocs = [IntLocs,locs2(ii):locs2(ii+1)-1];
        c = ii+1;
    end
end
IntLocs = IntLocs+1;
surf2_corr(IntLocs) = NaN;
surf2_corr(601:end) = interp1(601:length(surf2_corr),surf2_corr(601:end),601:length(surf2_corr),'pchip','extrap');

% Step 3: biggest interval of spikes
[~,locs3] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);
DiffSurf = diff(surf2_corr);
c = 0;
for ii = 1:length(locs3)-1
    if ii == c
        continue
    end
    if DiffSurf(locs3(ii))*DiffSurf(locs3(ii+1))<0 && abs(locs3(ii)-locs3(ii+1))<400
        IntLocs = [IntLocs,locs3(ii):locs3(ii+1)-1];
        c = ii+1;
    end
end
IntLocs = IntLocs+1;
surf2_corr(IntLocs) = NaN;
surf2_corr(601:end) = interp1(601:length(surf2_corr),surf2_corr(601:end),601:length(surf2_corr),'pchip','extrap');

%% Final single spikes
if strcmp(Type,'LFV')
    [~,locs4] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);locs4(locs4<1111 | locs4>9800) = [];
    locs4(locs4<1500) = [];
    locs4(locs4>3700) = [];
    for iiiii =1:length(locs4)
        surf2_corr(locs4(iiiii)-100:locs4(iiiii)+100) = NaN;
        surf2_corr(601:3700) = interp1(601:3700,surf2_corr(601:3700),601:3700,'spline','extrap');
    end
end

if strcmp(Type,'PIVSURF')
    [~,locs4] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);locs4(locs4<1111 | locs4>9800) = [];
    locs4(locs4<1500) = [];
    locs4(locs4>9000) = [];
    for iiiii = 1:length(locs4)
        Length = (locs4(iiiii)-500:locs4(iiiii)+500);
        SS = surf2_corr(Length);
        fit1 = polyfit(Length,SS,50);
        
        L_trans = 50;
        PolySurf2 = polyval(fit1,Length);
        MarkInt2(1) = locs4(iiiii)-25;
        MarkInt2(2) = locs4(iiiii)+25;
        X0 = locs4(iiiii)-(500+1);
        PS_mod = PolySurf2([MarkInt2(1):MarkInt2(end)]-X0);
        a = surf2_corr(1:MarkInt2(1)-(L_trans+1));
        c = surf2_corr(MarkInt2(end)+(L_trans+1):end);
        IndL = MarkInt2(1)-L_trans:MarkInt2(1);
        w1L = (MarkInt2(1)-[IndL])/(MarkInt2(1)-(MarkInt2(1)-L_trans));
        w2L = (IndL-(MarkInt2(1)-L_trans))/(MarkInt2(1)-(MarkInt2(1)-L_trans));
        bL = w1L.*surf2_corr(IndL)+w2L.*PolySurf2(IndL-X0);
        bL(end)= [];
        IndR = MarkInt2(end):MarkInt2(end)+L_trans;
        w1R = (IndR-(MarkInt2(end)))/(MarkInt2(end)+L_trans-MarkInt2(end));
        w2R = (MarkInt2(end)+L_trans-IndR)/(MarkInt2(end)+L_trans-MarkInt2(end));
        bR = w1R.*surf2_corr(IndR)+w2R.*PolySurf2(IndR-X0);
        bR(1) = [];

        surf2_corr = [a bL PS_mod bR c];
    end
end

end