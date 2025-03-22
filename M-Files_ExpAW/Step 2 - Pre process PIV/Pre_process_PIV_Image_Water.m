function [PIV] = Pre_process_PIV_Image_Water(PIV1)

PIV2 = adapthisteq(PIV1/max(PIV1(:)),'NBins',150)*1023;

H = 64; %32
B = sqrt(hanning(H).*hanning(H)');
B = B/sum(sum(B));
PIV = uint16(PIV2-filter2(B,PIV2,'same'));

PIV = double(PIV);

end
