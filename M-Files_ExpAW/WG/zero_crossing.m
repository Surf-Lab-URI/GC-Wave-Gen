function [Hraw,Z] = zero_crossing(x,tini,tfin,dt,noise_lev)

%       function H = zercro(x,dt,tini,tfin,noise_lev)
%
%       ---------------------------------------------
%       Zero-(up)crossing analysis of a wave train
%       ---------------------------------------------
%
%       input:
%             x        free surface elevation vector (can be a matrix of column vectors)
%             tini     [s]     starting time.
%             tfin	[s]	ending time
%             dt        [s]     intervallo tra due campioni
%             noise_lev         noise level (minimum height to consider)
%
%       output:  H              wave height vector (not ordered)
%                Z              vector = [nwaves,,Hrms,Hmed,Tmed,H3,T3,H10,T10,Hmax,Tmax] for each free surface column vector
%       output on screen:
%                 +-             -+
%                 | Hrms  no.onde |
%                 | Hmed  Tmed    |
%                 | H3    T3      |
%                 | H10   T10     |
%                 | Hmax  Tmax    |
%                 +-             -+
%

col=size(x,2);

for ii=1:col
    
    ita=x(:,ii);
    ntfin = fix(tfin / dt);
    if ntfin > length(ita)
        ntfin = length(ita);
    end
    
    ntini= fix(tini / dt);
    if ntini <= 0
        ntini = 1;
    end
    
    ita   = ita - mean(ita(ntini:1:ntfin));
    
    iin = ntini+1;
    while (ita(iin+1)-ita(iin)) <= 0 || (ita(iin)*ita(iin+1)) > 0
        iin = iin + 1;
        if iin > ntfin
            error('Not enough data')
        end
    end
    
    nonde = 0;
    tin   = dt / (ita(iin + 1) - ita(iin)) * ita(iin + 1);
    tcont = tin;
    ycrest= 0;
    ycav  = 0;
    
    %               loop sul numero di dati npun
    %               ----------------------------
    
    for i = iin:1:ntfin-1
        
        
        if ita(i) <= 0 && ita(i + 1) > 0 && (ycrest-ycav) > abs(noise_lev)
            nonde = nonde + 1;
            asup(nonde) = ycrest;
            ainf(nonde)= -ycav;
            H(nonde)=ycrest-ycav;
            tin = dt / (ita(i + 1) - ita(i)) * ita(i + 1);
            T(nonde) = dt - tin + tcont;
            tcont    = tin;
            ycrest   = 0.0;
            ycav     = 0.0;
        else
            tcont = tcont + dt;
            if ita(i) > ycrest
                ycrest = ita(i);
            end
            if ita(i) < ycav
                ycav = ita(i);
            end
        end
    end
    
    %               inizio conteggi zero-crossing
    %               sulle H e sui T complessivi
    
    Hrms = 0;
    for i = 1:nonde
        Hrms = Hrms + H(i) ^ 2;
    end
    Hrms = sqrt(Hrms / nonde);
    asuprms=sqrt(mean((asup.^2)));
    ainfrms=sqrt(mean((ainf.^2)));
    
    %
    % vettore altezza d'onda non ordinato
    %
    Hraw=H;
    
    %               Ordinamento in senso crescente
    
    [H,ind] = sort(H);
    T = T(ind(1:nonde));
    [asup,~] = sort(asup);
    [ainf,~] = sort(ainf);
    %               Ordinamento in senso decrescente
    
    H = H(nonde:-1:1);
    T = T(nonde:-1:1);
    asup=fliplr(asup);
    ainf=fliplr(ainf);
    
    %		Hmax, Tmax
    
    Hmax = H(1);
    Tmax = T(1);
    asupmax=asup(1);
    ainfmax=ainf(1);
    %               Tmed, Hmed
    
    Hmed = sum(H)/nonde;
    Tmed = sum(T)/nonde;
    asupmed=mean(asup);
    ainfmed=mean(ainf);
    %               H1/3, T1/3
    
    i3 = fix(nonde / 3);
    if i3 >= 1
        H3 = sum(H(1:i3))/i3;
        T3 = sum(T(1:i3))/i3;
        asup3=mean(asup(1:i3));
        ainf3=mean(ainf(1:i3));
    end
    
    %               H1/10, T1/10
    
    i10 = fix(nonde / 10);
    if i10 >= 1
        H10 = sum(H(1:i10))/i10;
        T10 = sum(T(1:i10))/i10;
        asup10=mean(asup(1:i10));
        ainf10=mean(ainf(1:i10));
        
    end
    
    %               H1/20, T1/20
    
    i20 = fix(nonde / 20);
    if i20 >= 1
        H20 = sum(H(1:i20))/i20;
        T20 = sum(T(1:i20))/i20;
    end
    
    %     y1(1)=Hrms;  y1(2)=Hmed; y1(3)=H3; y1(4)=H10; y1(5)=Hmax;
    % 	y2(1)=nonde; y2(2)=Tmed; y2(3)=T3; y2(4)=T10; y2(5)=Tmax;
    %
    %         Y = [y1(:) y2(:)];
    %         aasup=[ asuprms asupmed asup3 asup10 asupmax]';
    %         aainf=[ ainfrms ainfmed ainf3 ainf10 ainfmax]';
    Z(:,ii)=[nonde,Hrms,Hmed,Tmed,H3,T3,H10,T10,Hmax,Tmax];
end

