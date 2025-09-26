function [vector] = SUBPIXGAUSS (result_conv,interrogationarea,x,y,SubPixOffset)
if (x <= (size(result_conv,1)-1)) && (y <= (size(result_conv,1)-1)) && (x >= 1) && (y >= 1)
    %the following 8 lines are copyright (c) 1998, Uri Shavit, Roi Gurka, Alex Liberzon, Technion – Israel Institute of Technology
    %http://urapiv.wordpress.com
    f0 = log(result_conv(y,x));
    f1 = log(result_conv(y-1,x));
    f2 = log(result_conv(y+1,x));
    peaky = y+ (f1-f2)/(2*f1-4*f0+2*f2);
    f0 = log(result_conv(y,x));
    f1 = log(result_conv(y,x-1));
    f2 = log(result_conv(y,x+1));
    peakx = x+ (f1-f2)/(2*f1-4*f0+2*f2);
    %
    SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
    SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
    vector=[SubpixelX, SubpixelY];  
 else
    vector=[NaN NaN];
end