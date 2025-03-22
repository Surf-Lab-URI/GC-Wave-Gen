function surf2_corr = despike_jumps(surface2,Jump_thr,Type,expName)

%%%%%% Function used to despike Shoaling waves from sudden jumps in the
%%%%%% surface detection

% Step 0
surf2_corr = surface2;
L_poly = 1000;
L_trans = 150;
if strcmp(Type,'PIVSURF')
    Imin = 1500;
    Imax = 9300;
elseif strcmp(Type,'LFV')
    Imin = 801;
    Imax = 3900;
end

% Step 1: from small to large intervals of spikes
Int_Lev = [50 100 200 400 800];
for iii = 1:length(Int_Lev)
    [~,locs1] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);locs1(locs1<Imin | locs1>Imax) = [];
    DiffSurf = diff(surf2_corr);
    IntLocs = [];
    c = 0;
    for ii = 1:length(locs1)-1
        if ii == c 
            continue
        end
        if DiffSurf(locs1(ii))*DiffSurf(locs1(ii+1))<0 && abs(locs1(ii)-locs1(ii+1))<Int_Lev(iii)
            I1 = locs1(ii):locs1(ii+1)-1;
            Bound = [I1(1)-1 I1(end)+1];
            surf2_corr = vec_int(surf2_corr,Bound,L_poly,L_trans);
            IntLocs = [IntLocs,I1];
            c = ii+1;
        end
        IntLocs = IntLocs+1;
    end
end

%% Final single spikes
if strcmp(Type,'LFV')
    [~,locs4] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);locs4(locs4<1111 | locs4>9800) = [];
    
    if str2double(expName)>36 && str2double(expName)<49
        Lmin = 1100;
        Lmax = 3900;
    else
        Lmin = 1500;
        Lmax = 3700;
    end
    
    locs4(locs4<Lmin) = [];
    locs4(locs4>Lmax) = [];
    
    for iiiii =1:length(locs4)
        surf2_corr(locs4(iiiii)-100:locs4(iiiii)+100) = NaN;
        surf2_corr(Lmin-100:Lmax+100) = interp1(Lmin-100:Lmax+100,surf2_corr(Lmin-100:Lmax+100),Lmin-100:Lmax+100,'spline','extrap');
        %         surf2_corr(601:3700) = interp1(601:3700,surf2_corr(601:3700),601:3700,'spline','extrap');
    end
end

if strcmp(Type,'PIVSURF')
    [~,locs4] = findpeaks(abs(diff(surf2_corr)), 'minpeakheight', Jump_thr);locs4(locs4<1111 | locs4>9800) = [];
    locs4(locs4<1500) = [];
    locs4(locs4>9300) = [];
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