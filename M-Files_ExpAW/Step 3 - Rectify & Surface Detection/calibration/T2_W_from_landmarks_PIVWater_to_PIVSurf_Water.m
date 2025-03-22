%%% Find landmarks of PIV Water in PIVSurf Water coordinates 
% Points used for PIV Water coordinates references
%%% Projection of the grid points under the water
% Left points LP1 and LP2
LP_Inch8th = [1203 1963 ; 1203 1915 ; 1204 1866 ; 1204 1819 ; 1204 1772 ; 1205 1724 ; 1205 1675 ; 1206 1627 ; 1206 1579 ; 1207 1531 ; 1207 1483 ; 1207 1436 ; 1208 1387 ; 1208 1338 ; 1208 1290 ];
DeltaY_LP = -mean(diff(LP_Inch8th(:,2)));
CC1 = polyfit(LP_Inch8th(:,2),LP_Inch8th(:,1),1); % straight line fitting landmarks
LP1(1,2) = LP_Inch8th(1,2)+4*DeltaY_LP;
LP1(1,1) = polyval(CC1,LP1(1,2));
LP2(1,2) = LP1(1,2)+(3*8)*DeltaY_LP;
LP2(1,1) = polyval(CC1,LP2(1,2));
% Right points RP1 and RP2
RP_Inch8th = [ 3163 1986 ; 3163 1938 ; 3164 1890 ; 3164 1841 ; 3164 1793 ; 3165 1745 ; 3165 1697 ; 3165 1649 ; 3166 1601 ; 3166 1553 ; 3166 1503 ; 3167 1455 ;  3167 1407 ; 3168 1359 ;  3168 1311];
DeltaY_RP = -mean(diff(RP_Inch8th(:,2)));
CC2 = polyfit(RP_Inch8th(:,2),RP_Inch8th(:,1),1); % straight line fitting landmarks
RP1(1,2) = RP_Inch8th(1,2)+4*DeltaY_RP;
RP1(1,1) = polyval(CC2,RP1(1,2));
RP2(1,2) = RP1(1,2)+(3*8)*DeltaY_RP;
RP2(1,1) = polyval(CC2,RP2(1,2));


U2 = [ LP1 ; LP2 ; RP2 ; RP1]; % The coordinates of 
% four corner of a quadrilateral in the inbound image, or in the image that
% should be transformed.
X2 = [ 630 875 ; 610 2715 ; 3706 2743; 3718 903 ]; % The coordinates of 
% four corner of quadrilateral in the outbound image, or in the transformed
% image.

T2 = maketform('projective',U2,X2);