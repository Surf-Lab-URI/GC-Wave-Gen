% Script for generating images of the free surface to be assembled into a
% video. Use ffmpeg to make a video with a command like: ffmpeg -start_number 100  -framerate 3 -i ExpLCL_1_03_Pivsurf_%03d.jpg -frames:v 60 -vf "scale=1920:1080" -c:v libx264 -r 3 -pix_fmt yuv420p output2.mp4

clear
LONG = '/media/surflab/LC_Working24/LC/FabMarcNovDec2014/data/Longitudinal/PIVdt10ms_IRlas1_8hz/';
DIRS=dir(LONG);
DIRS=DIRS(3:end);


for ii=1%:length(DIRS)

    exp_name=DIRS(ii).name;
    
    num_of_digits = 3;
    load_path = [LONG exp_name];
    files=dir([load_path '/PIVRaw/PIV/*.mat']);
    number_of_pair=length(files)/2;
    
    
    image_pair_number = 100;
    previewing = true;
    for image_pair_number = 0:number_of_pair-1
    
        tic
        %PIV Surf
        im = load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_a.mat']); %replace ~ with path
        imgPivsurfa = im.imgPivsurf;%gpuArray(imgPivsurf); GPU array doesn't seem to make it any faster
        im = load([load_path '/PIVRaw/PIVSURF/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '_b.mat']); %replace ~ with path
        imgPivsurfb = im.imgPivsurf;%gpuArray(imgPivsurf); GPU array doesn't seem to make it any faster

        %Surface detection and Creating Masks
        U1 = [147 49;2024 57; 1995 1004; 161 999];
        X1 = [147 49; 2024 49; 2024 1004; 147 1004];
        T1 = fitgeotrans(U1,X1,'projective');
        da=imwarp(imgPivsurfa,T1,'cubic');
        db=imwarp(imgPivsurfb,T1,'cubic');
        
        da=imresize(da,176.9769/105.5880);%Resizing to match PIV
        db=imresize(db,176.9769/105.5880);%Resizing to match PIV
        ca = da(30:3525,150:end-150); %Cropping for surface video (Larger than PIV images)
        cb = db(30:3525,150:end-150); %Cropping for surface video (Larger than PIV images)
        sa=da(30:3525,755:755+2047); %cropping
        sb=db(30:3525,755:755+2047);
        
        %%%%%OG Fabrice surface detection
        % imSurf1 = findSurface_simple_ext_force_2023((medfilt2(s1)), 1);
        % imSurf2 = findSurface_simple_ext_force_2023((medfilt2(s2)), 1);
        
        %%%%%Andy Crapper-optimized surface detection
        surfSigmas = [50 40 30 20 15];
        surfSteps = [50 40 30 5];
        SurfMask = 1;
        slopeDiffThreshold = 5;
        imSurfa = CrapperOptimized_FindSurface(ca, surfSigmas, surfSteps, 1, slopeDiffThreshold);
        imSurfa.surface = FiltSurf(imSurfa.surface_raw,200);
        imSurfb = CrapperOptimized_FindSurface(cb, surfSigmas, surfSteps, 1, slopeDiffThreshold);
        imSurfb.surface = FiltSurf(imSurfb.surface_raw,200);
        
        
        disp(['image_pair_number = ', num2str(image_pair_number)])
        figure(1)
        hold off
        imagesc(ca, [130,300])
        hold on
        daspect([1,1,1])
        colormap gray
        axis off
        plot(imSurfa.surface,'-r','LineWidth',2)
        set(gca,'YLim', [1500,2100],'FontSize',20)
        set(gcf,'Units', 'Normalized', 'OuterPosition', [0.65078125,0,0.97421875,0.90625])

        % figure(2)
        % hold off
        % imagesc(cb, [130,300])
        % hold on
        % colormap gray
        % axis off
        % daspect([1,1,1])
        % plot(imSurfb.surface)
        % set(gca,'YLim', [1500,2100],'FontSize',20)
        % toc
        
        %Scale Bar
        figure(1) 
        xl = xlim;
        yl = ylim;
        
        lsbm = 1e-2; % length of scale bar in meters
        mpp = 1/17697.69; % meters per pixel
        lsb = lsbm/mpp;
        xsb = [xl(1)+(xl(2)-xl(1))*0.05, xl(1)+(xl(2)-xl(1))*0.05+lsb];
        ysb = (yl(2) - (yl(2)-yl(1))*0.1)*[1 1];
        try
            delete(sb)
            delete(sbt)
        end
        sb = plot(xsb,ysb,'-k', 'LineWidth',10);
    
    
        sbl = sprintf('%d cm',lsbm*100);
        sbt = text(xsb(2) + (xl(2)-xl(1))*0.01,ysb(2), sbl,'FontSize',24,'Interpreter','latex');

        text(0.5*xsb(1),(yl(2) - (yl(2)-yl(1))*0.35),'Water','FontSize',24,'Interpreter','latex')
        text(0.5*xsb(1),(yl(2) - (yl(2)-yl(1))*0.6),'Air','FontSize',24,'Interpreter','latex','Color',[1,1,1])
        text(0.5*xsb(1),(yl(2) - (yl(2)-yl(1))*0.75),'$Wind \rightarrow$','FontSize',24,'Interpreter','latex','Color',[1,1,1])

        ttext = sprintf('t = %.1f s',image_pair_number/7.2);
        text(xl(1)+(xl(2)-xl(1))*0.92, yl(2) - (yl(2)-yl(1))*0.94,ttext,'FontSize',24,'Interpreter','latex','Color',[1,1,1])
        drawnow
        
        % ip = input('a for back, d for forward','s');
        % nip = str2double(ip);
        % if ip == 'a'
        %     image_pair_number = max(0,image_pair_number-1);
        % elseif ip == 'd'
        %     image_pair_number = min(number_of_pair-1, image_pair_number+1);
        % elseif ~isnan(nip) && floor(nip) == nip && image_pair_number >= 0 && image_pair_number < number_of_pair
        %     image_pair_number = nip;
        % end

        fname = ['videoframes/' exp_name '_Pivsurf_' sprintf(['%0' num2str(num_of_digits) 'd'], image_pair_number) '.jpg']; % full name of image
        exportgraphics(gca,fname,'Resolution','1200')

    end
end

