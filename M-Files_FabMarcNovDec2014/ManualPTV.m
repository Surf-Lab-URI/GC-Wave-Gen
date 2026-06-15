clear
clc
% close all

LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
DIRS=dir(LONG);
DIRS=DIRS(3:end);

ii = 1;

exp_name=DIRS(ii).name;

num_of_digits = 3;
load_path = [LONG exp_name];
files=dir([load_path '/PIVRaw/PIV/*.mat']);
number_of_pair=length(files)/2;

image_pair_number = 123

load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
IM_a = imgPiv;
load([load_path '/PIVRaw/PIV/' exp_name '_Piv_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
IM_b = imgPiv;

imSurfa = FindSurfaceCapillary([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat'], findMask = true); 

imSurfb = FindSurfaceCapillary([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat'], findMask = true);



clicking = true;
np = 1;
i = 1;
p = nan(1,2,2);
ab = 1;
altmask_offset = 10;

% xl = [1, size(IM_a,2)];
% yl = [1, size(IM_a,1)];

xc = [300,800];
yc = [340,840];
IM_a_mask = IM_a .* imSurfa.mask;
IM_b_mask = IM_b .* imSurfb.mask;
IM_a = IM_a(yc(1):yc(2), xc(1):xc(2));
IM_b = IM_b(yc(1):yc(2), xc(1):xc(2));

imSurfa.MPTVSurf = imSurfa.surfacePIVImg(xc(1):xc(2))-yc(1);
imSurfb.MPTVSurf = imSurfb.surfacePIVImg(xc(1):xc(2))-yc(1);



xli = [1, size(IM_a,2)];
yli = [1, size(IM_a,1)];
xl = xli;
yl = yli;

while clicking

    figure(1)
    hold off
    if ab == 1
        imagesc(IM_a,[0,200])
        hold on
        plot(imSurfa.MPTVSurf,'-r')
        plot(imSurfa.MPTVSurf + altmask_offset,'-b')
        % xlim(xl)
        % ylim(yl)
        title(sprintf('Particle Num. %d,    Img. A', i))
        disp('select particle in image A')
    else
        imagesc(IM_b,[0,200])
        hold on
        plot(imSurfb.MPTVSurf,'-r')
        plot(imSurfb.MPTVSurf + altmask_offset,'-b')
        % colormap gray
        % daspect([1,1,1])
        % xlim(xl)
        % ylim(yl)
        title(sprintf('Particle Num. %d,    Img. B', i))
        disp('select particle in image B')
    end
    hold on
    daspect([1,1,1])
    colormap gray
    if i > 1
        plot(squeeze(p(1:end,ab,1)), squeeze(p(1:end,ab,2)), '+r','MarkerSize',5)
    end
    if size(p,1) >= i
        plot(squeeze(p(i,ab,1)), squeeze(p(i,ab,2)), '+g','MarkerSize',5)
    end

    
    [x,y,ip] = ginput(1);
    disp(ip)
    % nip = str2double(ip);
    if ip == 'a'
        ab = 1;
    elseif ip == 'd'
        ab = 2;
    elseif ip == 1
        p(i,ab,:) = [x,y];
    elseif ip == 'n'
        ip = input('a for back, d for forward','s');
        nip = str2double(ip);
        if ~isnan(nip) && floor(nip) == nip && image_pair_number >= 0 && image_pair_number < number_of_pair
            i = nip;
        end
    elseif ip == 's'
        if any(isnan(p(i,:,:)))
            disp('you have not selected a particle in both frames')
        else
            np = np + 1;
            i = np;
            ab = 1;
        end
    elseif ip == 'e'
        ipf = input('are you sure you are finished? (y,n)','s');
        if ipf == 'y'
            clicking = flase;
        end
    elseif ip == 'z'
        disp('select x extents')
        [x2,~] = ginput(2);
        disp('select y extents')
        [~,y2] = ginput(2);
        xl(1) = min(x2);
        xl(2) = max(x2);
        yl(1) = min(y2);
        yl(2) = max(y2);
    elseif ip =='o'
        xl = xli;
        yl = yli;
    elseif ip == 'x'
        altmask_offset = altmask_offset + 1;
    elseif ip == 'w'
        altmask_offset = altmask_offset - 1;
    end
end

if any(isnan(np,:,:))
    np = np-1;
    p = p(1:end-1,:,:);
end

imSurfa.altmask = zeros(size(IM_a));
imSurfb.altmask = zeros(size(IM_b));

for k = 1:length(imSurfa.MPTVSurf)
    imSurfa.altmask(imSurfa.MPTVSurf + altmask_offset:end, k) = 1;
    imSurfb.altmask(imSurfb.MPTVSurf + altmask_offset:end, k) = 1;
end
IM_a_altmask = IM_a .* imSurfa.altmask;
IM_b_altmask = IM_b .* imSurfb.altmask;

%% Write results
tifPath = 'GTImgs/';
if ~exist(tifPath, 'dir')
    mkdir(tifPath);
end
imwrite(uint8(IM_a),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgA.tif'])
imwrite(uint8(IM_b),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgB.tif'])
imwrite(uint8(IM_a_mask),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgA_mask.tif'])
imwrite(uint8(IM_b_mask),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgB_mask.tif'])
imwrite(uint8(IM_a_mask),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgA_altmask.tif'])
imwrite(uint8(IM_b_mask),[tifPath, exp_name,'-',int2str(image_pair_number), '_imgB_altmask.tif'])
surfa = imSurfa.MPTVSurf;
surfb = imSurfb.MPTVSurf;
outname = [tifPath, exp_name,'-',int2str(image_pair_number),'.mat'];
save(outname, "p","np","surfa", "surfb", "altmask_offset")


% if isfile(outfile)
%     ip = input("An output file for this pair already exists. Do you want to overwrite it? (Y/n)")
%     if ip ~= "Y"
%         ip2 = input("")
%     end
% end