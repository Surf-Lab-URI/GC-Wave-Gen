function matFrame = saveToMatSingleFrame(rawFrame,n,m)


fid = fopen(rawFrame);

img = fread(fid,n*m,'uint16');
img = reshape(img,n,m);
img = img';
tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 zeros
fclose(fid);
matFrame.ts = ts;
matFrame.img = img;
