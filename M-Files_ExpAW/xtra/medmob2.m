% function y = medmob(x,N)
%
% Media mobile di un segnale anche se contiene NaN
%
%               1  i=j+N/2 
%       y(j) =  -  Sommat.  x(i) 
%               N  i=j-N/2
%
% Input:
%	x	vettore dati (reali)
%	N	# di punti corrispondenti a T, su cui mediare
%
% Output:
%	y  	vettore dati (reali) calcolati
%
% N.B.		equivale ad un filtraggio con fn = 1 / (2*N*dt)
% 		ovvero: N = 1 / (2*fn*dt)
 
% By Marco Petti, DIC - FI.

% Modificata da F.Addona (08/02/2018): x può essere una matrice (vettori incolonnati)

function y = medmob2(x,N)

col = size(x,2);

for i = 1:col
    xx=x(:,i);
    NPUN = length(xx);
    N2   = fix(N/2);

    for j=1:N2
    %   y(j)=nanmean(x(1:j+N2));
       y(j,i)=nanmean(xx(1:j));

    end
    for j=NPUN-N2+1:NPUN
    %   y(j)=nanmean(x(NPUN-j+1:NPUN));
          y(j,i)=nanmean(xx(j:NPUN));

    end

    for j = N2+1:NPUN-N2,
        y(j,i) = nanmean(xx(j-N2:j+N2));
    end
    %y = y/(2*N2+1);

    %y(NPUN-N2+1:NPUN) = x(NPUN-N2+1:NPUN);
    %y(1:N2) = x(1:N2);

    %end
end
