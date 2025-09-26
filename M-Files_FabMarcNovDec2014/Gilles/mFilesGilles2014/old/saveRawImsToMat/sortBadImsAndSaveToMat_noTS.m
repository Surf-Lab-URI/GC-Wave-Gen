clear all
% close all
tic,
Dt = 1/7.2;
exp_name = 'LC1_dt25ms_1';
exp_folder = ['F:\data\Exp' exp_name '\rawImages'];
prefix = 'Movie10_Scene3_';
num_of_digits_lfv = 4;
num_of_digits_piv = 4;
% nC = 5; % number of 180 sec chunks
cL = 4319; %chunk length

%% Bad frame increments
bfiLfv = 0;
bfiPivsurf = 0;
bfiPiv1 = 0;
bfiPiv2 = 0;

%% for loop on chunks
cNum = 0;%:nC-1
    %% for loop on image number
    for imNum = 0:cL
%         imNum = 10;
        imNumR = imNum + cNum * (cL +1)
        lfvNum = floor(imNumR/2) + bfiLfv;
        pivsurfNum = imNumR + bfiPivsurf;
        piv1Num = imNumR + bfiPiv1;
        piv2Num = imNumR + bfiPiv2;
        %
        imLfv = [exp_folder '\Lfv\' prefix 'Lfv_' ...
            sprintf(['%0' num2str(num_of_digits_lfv) 'd'],lfvNum) '.raw'];
        imPivsurf = [exp_folder '\Pivsurf\' prefix 'Pivsurf_' ...
            sprintf(['%0' num2str(num_of_digits_piv) 'd'],pivsurfNum) '.raw'];
        imPiv1 = [exp_folder '\Piv1\' prefix 'Piv1_' ...
            sprintf(['%0' num2str(num_of_digits_piv) 'd'],piv1Num) '.raw'];
        imPiv2 = [exp_folder '\Piv2\' prefix 'Piv2_' ...
            sprintf(['%0' num2str(num_of_digits_piv) 'd'],piv2Num) '.raw'];
        n = 2048; %number of rows
        m = 2048; %number of columns
        
        %% Lfv
        if mod(imNumR,2)==0
            fid = fopen(imLfv);
            %
            img = fread(fid,n*m,'uint16');
            img = reshape(img,n,m);
            img = img';
%             tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%             ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 zeros
            fclose(fid);
            %
            imgLfv = rot90(img,2);  %lfv cam is upside and flipped left-right
%             tsLfv = ts;
        else
            imgLfv = nan;
%             tsLfv = nan;
        end
        %% Pivsurf
        fid = fopen(imPivsurf);
        %
        img = fread(fid,n*m,'uint16');
        img = reshape(img,n,m);
        img = img';
%         tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%         ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
        fclose(fid);
        %
        imgPivsurf = flipud(img);  %pivsurf camera is upside down
%         tsPivsurf = ts;
        %% Piv1
        fid = fopen(imPiv1);
        %
        img = fread(fid,n*m,'uint16');
        img = reshape(img,n,m);
        img = img';
%         tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%         ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
        fclose(fid);
        %
        imgPiv1 = fliplr(img); %Piv1 cam is flipped lr
%         tsPiv1 = ts;
        %% Piv2
        fid = fopen(imPiv2);
        %
        img = fread(fid,n*m,'uint16');
        img = reshape(img,n,m);
        img = img';
%         tsTemp1 = fread(fid, 'uint64'); %timestamp is a footer
%         ts = tsTemp1(1);  %tsTemp1 is a vector padded with numofFrames-1 seros
        fclose(fid);
        %
        imgPiv2 = fliplr(img); %Piv2 cam is flipped lr
%         tsPiv2 = ts;
        %% Tests
%         tss = [tsLfv;  tsPivsurf; tsPiv1; tsPiv2];
        nums = [lfvNum; pivsurfNum ; piv1Num ; piv2Num ];
%         if imNum == 0
% %             tss
%             nums
%             pause
% %             save(['ts0_exp' exp_name '.mat'], 'tss');
%         end
        %         if tsPivsurf~=tsPiv1 || tsPivsurf~=tsPiv2 || tsPiv1~=tsPiv2
        
        % ts = [tsPivsurf; tsPiv1; tsPiv2]
%         ts0 = load(['ts0_exp' exp_name '.mat']);
        % ts0_3 = ts0(2:end);
%         dif = imNum/2 - (tss-ts0.tss)*1d-7/Dt;
%         pause
%         if mod(imNum,2)==0 && sum(abs(dif)>3d-2)>0 || mod(imNum,2)==1 && sum(abs(dif)>5d-2)>0
%             dif
%             pause
%         end
        %% save name example: Exp3_Piv2_0440_a
        if mod(imNumR,2) == 0
            imgLfvSaveName = [exp_folder '\Lfv\Exp' exp_name '_Lfv_' ...
                sprintf(['%0' num2str(num_of_digits_lfv) 'd'],floor(imNumR/2)) '.mat'];
            eval(['save ' imgLfvSaveName ' imgLfv']);
            pivSaveNum = [sprintf(['%0' num2str(num_of_digits_piv) 'd'],floor(imNumR/2)) '_a'];
        else
            pivSaveNum = [sprintf(['%0' num2str(num_of_digits_piv) 'd'],floor(imNumR/2)) '_b'];
        end
        %
        imgPivsurfSaveName = [exp_folder '\Pivsurf\Exp' exp_name '_Pivsurf_' ...
            pivSaveNum '.mat'];
        eval(['save ' imgPivsurfSaveName ' imgPivsurf']);
        %
        imgPiv1SaveName = [exp_folder '\Piv1\Exp' exp_name '_Piv1_' ...
            pivSaveNum '.mat'];
        eval(['save ' imgPiv1SaveName ' imgPiv1']);
        %
        imgPiv2SaveName = [exp_folder '\Piv2\Exp' exp_name '_Piv2_' ...
            pivSaveNum '.mat'];
        eval(['save ' imgPiv2SaveName ' imgPiv2']);
        %
    end
toc

