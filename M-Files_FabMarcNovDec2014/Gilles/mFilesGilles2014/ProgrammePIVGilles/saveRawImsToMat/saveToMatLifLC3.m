clear all
close all

%% user input
exp_name = 'LC1_dt15ms_cc_LIF_1';
prefix = 'Movie15_Scene5';
%%
Dt = 1/7.2;
exp_folder = ['\\beo\data\Exp' exp_name '\'];
rawImFolder = ['\\beo\data\rawPIV\' prefix];
num_of_digits_lfv = 4;
num_of_digits_piv = 4;
nI = 2591; % number of images - 1

%% Bad frame increments
bfiLfv = 0;
bfiPivsurf = 0;
bfiPiv1 = 0;
bfiPiv2 = 0;

    %% for loop on image number
    for imNum = 0:nI
%         imNum = 0;
%         imNumR = imNum + cNum * (cL +1)
        lfvNum = floor(imNum/2) + bfiLfv;
        pivsurfNum = lfvNum;%imNum + bfiPivsurf;
        piv1Num = imNum + bfiPiv1;
        piv2Num = imNum + bfiPiv2;
        %
%         imLfv = [rawImFolder '\Lfv\' prefix '_lfv_' ...
%             sprintf(['%0' num2str(num_of_digits_lfv) 'd'],lfvNum) '.raw'];
%         imPivsurf = [rawImFolder '\Pivsurf\' prefix '_pivsurf_' ...
%             sprintf(['%0' num2str(num_of_digits_piv) 'd'],pivsurfNum) '.raw'];
%         imPiv1 = [rawImFolder '\Pivsurfcc\' prefix '_pivsurfcc_' ...
%             sprintf(['%0' num2str(num_of_digits_piv) 'd'],piv1Num) '.raw'];
        imPiv2 = [rawImFolder '\Pivcc\' prefix '_pivcc_' ...
            sprintf(['%0' num2str(num_of_digits_piv) 'd'],piv2Num) '.raw'];
        n = 2048; %number of rows
        m = 2048; %number of columns
        
%         %% Lfv
%         if mod(imNum,2)==0
%             fid = fopen(imLfv);
%             %
%             img = fread(fid,n*m,'uint16');
%             img = reshape(img,n,m);
%             img = img';
%             tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%             ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 zeros
%             fclose(fid);
%             %
%             imgLfv = rot90(img,2);  %lfv cam is upside and flipped left-right
%             tsLfv = ts;
%         else
%             imgLfv = nan;
%             tsLfv = nan;
%         end
%         %% Pivsurf
%         fid = fopen(imPivsurf);
%         %
%         img = fread(fid,n*m,'uint16');
%         img = reshape(img,n,m);
%         img = img';
%         tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%         ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
%         fclose(fid);
%         %
%         imgPivsurf = flipud(img);  %pivsurf camera is upside down
%         tsPivsurf = ts;
%         %% Piv1
%         fid = fopen(imPiv1);
%         %
%         img = fread(fid,n*m,'uint16');
%         img = reshape(img,n,m);
%         img = img';
%         tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%         ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
%         fclose(fid);
%         %
%         imgPivsurfcc = rot90(img,2); %Pivsurfcc
%         tsPivsurfcc = ts;
        %% Piv2
        fid = fopen(imPiv2);
        %
        img = fread(fid,n*m,'uint16');
        img = reshape(img,n,m);
        img = img';
        tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
        ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
        fclose(fid);
        %
        imgPivcc = img; %Pivcc
        tsPivcc = ts;
        %% Tests
%         tss = [tsLfv;  tsPivsurf; tsPivsurfcc; tsPivcc];
%         nums = [lfvNum; pivsurfNum ; piv1Num ; piv2Num ];
%         if imNum == 0
%             tss
%             nums
%             pause
%             save(['ts0_exp' exp_name '.mat'], 'tss');
%         end
        %         if tsPivsurf~=tsPivsurfcc || tsPivsurf~=tsPivcc || tsPivsurfcc~=tsPivcc
        
        % ts = [tsPivsurf; tsPivsurfcc; tsPivcc]
%         ts0 = load(['ts0_exp' exp_name '.mat']);
        % ts0_3 = ts0(2:end);
%         dif = imNum/2 - (tss-ts0.tss)*1d-7/Dt;
%         pause
%         if mod(imNum,2)==0 && sum(dif(:)>7.5d-2)>0 || mod(imNum,2)==1 && sum(dif(:)>7.5d-2)>0% || sum(dif<0)>0
%             dif
%             pause
%         end
        %% save name example: Exp3_Piv2_0440_a
        if mod(imNum,2) == 0
%             imgLfvSaveName = [exp_folder 'RawImages\Lfv\Exp' exp_name '_Lfv_' ...
%                 sprintf(['%0' num2str(num_of_digits_lfv) 'd'],floor(imNum/2)) '.mat'];
%             imgPivsurfSaveName = [exp_folder 'RawImages\Pivsurf\Exp' exp_name '_Pivsurf_' ...
%                 sprintf(['%0' num2str(num_of_digits_lfv) 'd'],floor(imNum/2)) '.mat'];
%             eval(['save ' imgLfvSaveName ' imgLfv tsLfv']);
%             eval(['save ' imgPivsurfSaveName ' imgPivsurf tsPivsurf']);
            pivSaveNum = [sprintf(['%0' num2str(num_of_digits_piv) 'd'],floor(imNum/2)) '_a'];
        else
            pivSaveNum = [sprintf(['%0' num2str(num_of_digits_piv) 'd'],floor(imNum/2)) '_b'];
        end
        %
%         imgPivsurfSaveName = [exp_folder 'RawImages\Pivsurf\Exp' exp_name '_Pivsurf_' ...
%             pivSaveNum '.mat'];
%         eval(['save ' imgPivsurfSaveName ' imgPivsurf tsPivsurf']);
        %
%         imgPivsurfccSaveName = [exp_folder 'RawImages\Pivsurfcc\Exp' exp_name '_Pivsurfcc_' ...
%             pivSaveNum '.mat'];
%         eval(['save ' imgPivsurfccSaveName ' imgPivsurfcc tsPivsurfcc']);
        %
        imgPivccSaveName = [exp_folder 'RawImages\Pivcc\Exp' exp_name '_Pivcc_' ...
            pivSaveNum '.mat'];
        eval(['save ' imgPivccSaveName ' imgPivcc tsPivcc']);
        %
    end


