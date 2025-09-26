% Function fuseImages
%% Object:
% concatenates 2 PIV images (PIV1 and PIV2) together, with slight
% overlap, with linear gradient at overlap region AND flips final fused
% image left-right, for wind to come from the left
%% Arguments: 
% experiment number, paths of data, image pair number, number of digits in image pair number, image letter
%% Result: 
% experiment number, image pair number, image letter, fused image
%% Author:
% Marc Buckley
%% Last update:
% 04/29/2013
%% Example:
% exp_name = '3';
% path = ['\\Afsx1\piv2\Exp' exp_name '\' 'test_ims\'];
% image_pair_number = 440;
% num_of_digits = 4;
% image_letter = 'a';
% fusedIm = fuseImages(exp_name, path, image_pair_number, num_of_digits, image_letter);
% ans = 
%        fused_im: [2040x3944 double]
%        exp_name: '3'
%     im_pair_num: 440
%       im_letter: 'a'
%%
function fusedIm = fuseImages(exp_name, path1, path2, image_pair_number, num_of_digits, image_letter)
%[fusedIm.exp_name, fusedIm.im_pair_num, fusedIm.im_letter, fusedIm.fused_im]

piv1_struc = load([path1 'Exp' exp_name '_Piv1_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
piv1 = fliplr(piv1_struc.imgPiv1);
%
piv2_struc = load([path2 'Exp' exp_name '_Piv2_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_' image_letter]);
piv2 = fliplr(piv2_struc.imgPiv2);
%
d1 = 2048;
ncc =152;  % number of commun columns
nvc = 9;  % number of vertical shift pixels

c = piv1(1:end-nvc+1,ncc+1:end); %lower fraction of piv1 is kept

a = piv2(nvc:end,1:d1-ncc); %upper fraction of piv2 is kept, piv2 is downwind (left) of piv1

[s1,~] = size(a);

%middle piece: fusion of right edge of piv2 and left edge of piv1
b = nan(s1,ncc);
cc1 = piv1(1:end-nvc+1,1:ncc);
cc2 = piv2(nvc:end,d1 - ncc + 1:end);
% b = (cc1+cc2)/2; %optional brute force average, instead of following 4
% lines

%fuse with linear gradient of each image vertical border strip
for i = 1:ncc
    b(:,i) = ((ncc-i) * cc2(:,i) + i * cc1(:,i)) ./ ncc;
end

% concatenate 3 pieces
fusedIm.fused_im = fliplr([a b c]);
% save other info
fusedIm.exp_name = exp_name;
fusedIm.im_pair_num = image_pair_number;
fusedIm.im_letter = image_letter;