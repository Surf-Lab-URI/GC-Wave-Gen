function [PIV] = Pre_process_PIV_Image_Air_OLD(PIV1)

PIV2 = PIV1;
PIV1(PIV1==0)=NaN;
PIV1 = PIV1-min(min(PIV1));
PIV1 = PIV1/max(max(PIV1))*1023;

H = 64; %32
B = sqrt(hanning(H).*hanning(H)');
B = B/sum(sum(B));
PIV = PIV1-filter2(B,PIV1,'same');

PIV = PIV-min(min(PIV));
PIV = PIV/max(max(PIV))*1023;

end
