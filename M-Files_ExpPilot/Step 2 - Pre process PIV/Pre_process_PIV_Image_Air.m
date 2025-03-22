function [PIV] = Pre_process_PIV_Image_Air(PIV1)

PIV1 = PIV1-min(min(PIV1)); %Offset by minimum value
PIV1 = PIV1/max(max(PIV1))*1023; % Scale to 10 bits based on max value

H = 64; %32
B = sqrt(hanning(H).*hanning(H)'); %Hann matrix for smoothing
B = B/sum(sum(B)); %Normalize by total
PIV = PIV1-filter2(B,PIV1,'same'); %Smooth the image using 64x64 pixel windows

%redo offset and scale
PIV = PIV-min(min(PIV));
PIV = PIV/max(max(PIV))*1023;

end
