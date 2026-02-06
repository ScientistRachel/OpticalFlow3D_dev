clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\AiryScan_twoChannels\';
saveDir = [imDir '\PublicationFigures\Figure3\'];
movieDir = [imDir '\PublicationFigures\MovieS3\'];
ch0Name = 'C2-04302024_U2OS_Cherry-Myo_Lifeact-GFP_30s_3slices_1-Airyscan Processing-07'; % myosin
ch1Name = 'C1-04302024_U2OS_Cherry-Myo_Lifeact-GFP_30s_3slices_1-Airyscan Processing-07'; % actin
ch0Dir = [imDir filesep 'OpticalFlow3D' filesep ch0Name];
ch1Dir = [imDir filesep 'OpticalFlow3D' filesep ch1Name];

% Metadata
xyscale = 0.0425177; % um/pixel
zscale = 0.3; % um/pixel
tscale = 30/60; % minutes/frame
tRange = 4:28;
crop = 0; % No cropping necessary for this movie

% reliability thresholding
relPer = 40;

%%% Plotting parameters
gap = 2;
qscale = 1;
roi = [329 377 317 365]; % in pixels [x1 x2 y1 y2]

if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
if ~exist(movieDir,'dir')
    mkdir(movieDir)
end

%% Get Image Sizes

relFile = [ch0Name '_rel_t' num2str(4,'%04u') '.tiff'];
meta = imfinfo([ch0Dir filesep relFile]);
Nx = meta(1).Width;
Nx = Nx-2*crop;
Ny = meta(2).Height;
Ny = Ny-2*crop;
Nz = length(meta);

x = (0:Nx-1)*xyscale;
y = (0:Ny-1)*xyscale;
z = (0:Nz-1)*zscale;
[X,Y,Z] = meshgrid(x,y,z);
clear x y z

xD = X(1:gap:end);
yD = Y(1:gap:end);
% zD = Z(1:gap:end);


%% Load in images

myo = TIFFvolume([imDir filesep ch0Name '.tif'],31*3);
myo = reshape(myo,[size(myo,1), size(myo,2), 3, 31]);
myo = myo(crop+1:end-crop,crop+1:end-crop,:,:);
A = squeeze(max(myo,[],3));

act = TIFFvolume([imDir filesep ch1Name '.tif'],31*3);
act = reshape(act,[size(act,1), size(act,2), 3, 31]);
act = act(crop+1:end-crop,crop+1:end-crop,:,:);
B = squeeze(max(act,[],3));


%% Loop through frames
for exampleFrame = tRange

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Scale bar and time stamp on images

    imA = A(:,:,exampleFrame);
    imA = imadjust(imA,[30 2500]/2^16);

    if exampleFrame == 4
        a = size(imA,1)-30;
        b = size(imA,2)-30;
        imA_Fig4A = insertShape(imA,'Line',[b-5/xyscale a b a],'Color',2^16*[1 1 1],'LineWidth',20); % scalebar
        imwrite(imA_Fig4A,[saveDir filesep 'Fig3A_Myosin_Scale5um.png'])
        imA_Fig4A = insertShape(imA_Fig4A,'Rectangle',[roi(1) roi(3) roi(2)-roi(1) roi(4)-roi(3)],'Color',[204, 121, 167]*(2^16-1)/255,'LineWidth',5); % roi
        imwrite(imA_Fig4A,[saveDir filesep 'Fig3A_Myosin_Scale5um_withROI.png'])
    end

    imA = uint8(double(imA)/(2^16-1)*255);

    imB = B(:,:,exampleFrame);
    imB = imadjust(imB,[30 2500]/2^16);
    if exampleFrame == 4
        imB_Fig4A = insertShape(imB,'Line',[b-5/xyscale a b a],'Color',2^16*[1 1 1],'LineWidth',20); % scalebar
        imwrite(imB_Fig4A,[saveDir filesep 'Fig3A_Actin_Scale5um.png'])
    end
    imB = uint8(double(imB)/(2^16-1)*255);
    imB = repmat(imB,[1 1 3]);

    % Add timestamp for movie
    tText = (exampleFrame-1)*tscale;
    mm = floor(tText);
    ss = tText-mm;
    ss = round(ss*60); % seconds
    tText = [num2str(mm,'%02u') ':' num2str(ss,'%02u')];
    imA = insertText(imA,[20 15],tText,BoxOpacity=0,TextColor='white',FontSize=60);


    %% Reliability Thresholding
    
    relFile = [ch0Name '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    rel = TIFFvolume([ch0Dir filesep relFile],Nz);
    rel = rel(crop+1:end-crop,crop+1:end-crop,:);
    % disp('rel loaded')

    relThresh = prctile(rel(:),relPer);
    % disp(['Reliability Threshold Ch0: ' num2str(relThresh)])

    relMask0 = rel > relThresh;    
   
    relFile1 = [ch1Name '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    rel1 = TIFFvolume([ch1Dir filesep relFile1],Nz);
    rel1 = rel1(crop+1:end-crop,crop+1:end-crop,:);
    % disp('rel loaded')
   
    relThresh1 = prctile(rel1(:),relPer);
    % disp(['Reliability Threshold Ch1: ' num2str(relThresh1)])

    relMask1 = rel1 > relThresh1;

    relMask = relMask0 | relMask1; % Combine them for a fair comparison of dot products, magnitude, etc.
    % Using the OR to be somewhat relaxed with this constratint.
    % Technically any rel > 0 is a valid solution, so set that as a bare
    % minimum.
    relMin0 = rel > 0;
    relMin1 = rel1 > 0;
    relMask = relMask & relMin0 & relMin1;
    relMask = bwareaopen(relMask,10^5); % only keep the cell sized object

    %%% Clean up
    clear relMask0 relMask1

    %% Load in velocities
    
    vx = TIFFvolume([ch0Dir filesep ch0Name '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx = vx(crop+1:end-crop,crop+1:end-crop,:);
    vx = vx.*relMask./relMask*xyscale/tscale;

    vy = TIFFvolume([ch0Dir filesep ch0Name '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy = vy(crop+1:end-crop,crop+1:end-crop,:);
    vy = vy.*relMask./relMask*xyscale/tscale;

    vx1 = TIFFvolume([ch1Dir filesep ch1Name '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx1 = vx1(crop+1:end-crop,crop+1:end-crop,:);
    vx1 = vx1.*relMask./relMask*xyscale/tscale;

    vy1 = TIFFvolume([ch1Dir filesep ch1Name '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy1 = vy1(crop+1:end-crop,crop+1:end-crop,:);
    vy1 = vy1.*relMask./relMask*xyscale/tscale;

    theta = atan2(vy,vx);   
    theta1 = atan2(vy1,vx1);

    %% Theta Quiver Ch0

    vxD = vx(1:gap:end);
    vyD = vy(1:gap:end);
    vD = theta(:,:,2);
    vD = vD(1:gap:end);

    bins = linspace(-pi-eps,pi+eps,15);
    cmap = colorcet('C8');
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    figure(3)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD>=bins(jj)) & (vD<bins(jj+1));
        j = quiver(xD(slice),yD(slice),vxD(slice),vyD(slice),0,...
            'Color',cmap(jj,:),'LineWidth',0.25);
        hold on
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    set(gcf,'Color','w')
    set(gca,'XTick',[],'YTick',[])
    xlim([0 (Nx-1)*xyscale])
    ylim([0 (Ny-1)*xyscale])
    drawnow;

    quivA = getframe(gca);

    if exampleFrame == 4
        imwrite(quivA.cdata,[saveDir filesep 'Fig3B_Myosin_ThetaQuiver.png'])
    end

    quivA = imresize(quivA.cdata,[size(imA,1) size(imA,2)]);

    %% Theta Quiver Ch0 - ROI

    if exampleFrame == 4 % This is not part of the movie, just a figure panel

        imROI = imA(roi(3):roi(4),roi(1):roi(2));
        figure(31)
        imshow(imROI)
        imwrite(imROI,[saveDir filesep 'Fig3C_Myosin_ROI.png'])


        vxD = vx(1:gap:end);
        vyD = vy(1:gap:end);
        vD = theta(:,:,2);
        vD = vD(1:gap:end);
    
        bins = linspace(-pi-eps,pi+eps,15);
        cmap = colorcet('C8');
        cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);
    
        figure(32)
        set(gcf,'Position',[300 150 1250 1020])
        for jj = 1:length(bins)-1     
            slice = (vD>=bins(jj)) & (vD<bins(jj+1));
            j = quiver(xD(slice),yD(slice),vxD(slice),vyD(slice),0,...
                'Color',cmap(jj,:),'LineWidth',2);
            hold on
            if qscale ~= 0
                hU = get(j,'UData') ;
                hV = get(j,'VData') ;
                set(j,'UData',qscale*hU,'VData',qscale*hV)
            end
        end
        hold off
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'ydir','reverse')
        set(gca,'Color','k')
        set(gcf,'Color','w')
        set(gca,'XTick',[],'YTick',[])
        % xlim([14 16])
        % ylim([13.5 15.5])
        xlim(roi(1:2)*xyscale)
        ylim(roi(3:4)*xyscale)
        drawnow;
    
        imwrite(getframe(gca).cdata,[saveDir filesep 'Fig3C_Myosin_ThetaQuiver_ROI.png'])

    end

    %% Theta Quiver Ch1

    vxD = vx1(1:gap:end);
    vyD = vy1(1:gap:end);
    vD = theta1(:,:,2);
    vD = vD(1:gap:end);

    figure(4)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD>=bins(jj)) & (vD<bins(jj+1));
        j = quiver(xD(slice),yD(slice),vxD(slice),vyD(slice),0,...
            'Color',cmap(jj,:),'LineWidth',0.25);
        hold on
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    % set(gcf,'Color','w')
    set(gca,'XTick',[],'YTick',[])
    xlim([0 (Nx-1)*xyscale])
    ylim([0 (Ny-1)*xyscale])
    drawnow;

    quivB = getframe(gca);

    if exampleFrame == 4
        imwrite(quivB.cdata,[saveDir filesep 'Fig3B_Actin_ThetaQuiver.png'])
    end

    quivB = imresize(quivB.cdata,[size(imA,1) size(imA,2)]);

    quivB = insertShape(quivB,'Line',[b-5/xyscale a b a],'Color',2^16*[1 1 1],'LineWidth',20); % scalebar

    %% Combine Movie Panels

    figure(5)
    imshow([imA quivA ; imB quivB])
    imwrite([imA quivA ; imB quivB],[movieDir filesep 'MovieS3_frame' num2str(exampleFrame,'%03u') '_scaleBar5um_clockMMSS.png'])
   

    %% Clean up
    % close all
    disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end

%% Clean up
disp('Movie Complete')
close all