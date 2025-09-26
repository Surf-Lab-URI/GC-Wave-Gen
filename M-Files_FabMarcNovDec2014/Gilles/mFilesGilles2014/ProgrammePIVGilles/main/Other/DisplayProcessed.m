clear all
close all
exp_name = 'LC2_dt7ms_3';
calc_name = '_MOTDec';
start_image_number = 400;
end_image_number = 400;
path = '\\beo\data\';
num_of_digits = 4;
save_path = ['\\beo\data\Exp' exp_name '\ProcessedVelocities\'];

for image_pair_number=start_image_number:end_image_number
    disp(image_pair_number);
    %     load(['\\beo\data\Exp' exp_name '\ComputedVelocities\Exp' exp_name '_compVel_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '']);
    load(['\\beo\data\Exp' exp_name '\ProcessedVelocities\Exp' exp_name calc_name sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number)]);
%     figure, imagesc(MOTDec.up), colorbar
%     figure, imagesc(MOTDec.vp), colorbar
    figure, plot(MOTDec.up2m_rs)
    hold on, plot(MOTDec.vp2m_rs,'r')
    hold on, plot(MOTDec.upvpm_rs,'g')
    hold on, plot(2*MOTDec.ecpm_rs, 'k')
end