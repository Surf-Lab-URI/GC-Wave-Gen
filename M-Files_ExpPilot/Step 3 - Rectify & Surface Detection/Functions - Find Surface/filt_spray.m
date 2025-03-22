function [s_int] = filt_spray(ps_surf_raw)


s = ps_surf_raw;
%  s1 = s;
ds = diff(s); % Differences and approximate derivatives
b = 1;
a = 1;

Step = 500;
PW = 20; % diff peak width
TH = 50; % diff peak treshold

while a+Step < length(s)
    
    %     s(diff(ps_surf)>50) = nan;
    ds_rest = ds(b:end);
    a0 = find(abs(ds_rest)>TH,1, 'first') ;
    
    %     a_neg = find(ds(B:end)<-100*mean(abs(ds)),1, 'first');
    
    if isempty(a0);
        break;
    end
    
    a1 = a0 + b - 1;
    segi = ds(a1+PW:end);
    
    if ds(a1) > 0
        b0 = find(segi < -TH, 1, 'first');
    else
        b0 = find(segi > +TH, 1, 'first');
    end
    
    if isempty(b0)&&a1<500 % bad stuff is at the edge of image
        a = 1;
        b = a1 + PW;
    elseif isempty(b0)&&length(s)-a1<500
        a = a1 - PW;
        b = length(s);
    elseif isempty(b0)||(b0-a0)>500
        break;
    else
        b1 = a1 + b0 -1 + PW;
        b = b1 + PW;
        a = a1 - PW;
        %     keyboard,
    end
    
     if b>length(ps_surf_raw)
         b=length(ps_surf_raw);
     end
     if a<1
         a=1;
     end
    s(a:b) = nan;
    
end

x1 = find(~isnan(s));
s_int = interp1(x1,s(x1), 1:length(s), 'linear');
s_int_extrap = interp1(x1,s(x1), 1:length(s), 'nearest', 'extrap');
s_int(isnan(s_int)) = s_int_extrap(isnan(s_int)); % Replace ends with nearest extrapolated 

end

%
% a = find(ds(b+a+10>100*mean(abs(ds)),1, 'first');
%
% ds = diff(s);
% ds1 = ds;
% ds1(abs(ds)>100*mean(abs(ds))) = nan;
% figure, plot(ds)
% figure, plot(s(2:end));
% s1 = nancumsum(ds1,1,1);
%
% hold on, plot(s1, 'r')