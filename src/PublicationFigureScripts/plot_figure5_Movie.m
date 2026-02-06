clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\SiMView_Drosophila';
dataName = 'fused_tp__ch_0';
dataDir = [imDir filesep 'OpticalFlow3D' filesep dataName];
% saveDir = 'X:\Force Project\PublicationData\SiMView_Drosophila\PublicationFigures\MovieS7';
% figDir = 'X:\Force Project\PublicationData\SiMView_Drosophila\PublicationFigures\Figure5';
saveDir7 = 'Y:\instruments\ForceProject\MovieS7';
saveDir8 = 'Y:\instruments\ForceProject\MovieS8';

% Metadata
xyscale = 0.4114940; % um/pixel
zscale = 1.4963417;
tscale = 30/60; % minutes/frame

% reliability thresholding
relPer = 85;
minSize = 10^5;

% For plots
qscale = 1; % vector scaling, 1 or 0 both equal actual displacement
gap = 2; % show every vector or only every gap'th vector?
apMax = 2.5; % um/min
dvMax = 1.25; % um/min

% ROI -- remove empty padding around the embryo
x1 = 210;
x2 = 818;
y1 = 250;
y2 = 1700;
z1 = 5;
z2 = 160;

% clean figures
cleanFig = true;

%% Set up the save directories

if ~exist(saveDir7,'dir')
    mkdir(saveDir7)
end
if ~exist(saveDir8,'dir')
    mkdir(saveDir8)
end

%% Set up general matrices rather than recalculating each time
relFile = [dataName '_rel_t' num2str(4,'%04u') '.tiff'];
meta = imfinfo([dataDir filesep relFile]);
Nx = meta(1).Width;
Ny = meta(2).Height;
Nz = length(meta);

x = (0:Nx-1)*xyscale;
y = (0:Ny-1)*xyscale;
z = (0:Nz-1)*zscale;
[X,Y,Z] = meshgrid(x,y,z);
clear x y z

X = X(y1:y2,x1:x2,z1:z2);
Y = Y(y1:y2,x1:x2,z1:z2);
Z = Z(y1:y2,x1:x2,z1:z2);

%% Loop through frames
for exampleFrame = 4:88

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Load Images

    im1_ch0 = TIFFvolume([imDir filesep 'fused_tp_' num2str(exampleFrame-1) '_ch_0.tif'],Nz);
    % disp('raw images loaded')

    %% Reliability Thresholding
    
    relFile = [dataName '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    rel = TIFFvolume([dataDir filesep relFile],Nz);
    rel = rel(y1:y2,x1:x2,z1:z2);
    % disp('rel loaded')

    relThresh = prctile(rel(:),relPer);
    % disp(['Reliability Threshold: ' num2str(relThresh)])
    
    relMask = rel > relThresh;
    relMask = bwareaopen(relMask,minSize);

    % imagesc(sum(relMask,3)), set(gca,'DataAspectRatio',[1 1 1])
    % for kk = 1:size(relMask,3)
    %     imshow(relMask(:,:,kk));
    %     drawnow;
    % end
    
    %% Load in velocities
    
    vx = TIFFvolume([dataDir filesep dataName '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx = vx(y1:y2,x1:x2,z1:z2);
    vx = vx.*relMask./relMask*xyscale/tscale;    
    % disp('vx loaded')
    
    vy = TIFFvolume([dataDir filesep dataName '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy = vy(y1:y2,x1:x2,z1:z2);
    vy = vy.*relMask./relMask*xyscale/tscale;
    % disp('vy loaded')
    
    vz = TIFFvolume([dataDir filesep dataName '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vz = vz(y1:y2,x1:x2,z1:z2);
    vz = vz.*relMask./relMask*zscale/tscale;
    % disp('vz loaded')
    
    %% Shared quiver values

    xD = X(1:gap:end,1:gap:end,:);
    xD = xD-min(xD(:)); % Zeroing here lets me match the quiver and image figure camera positions effectively.
    yD = Y(1:gap:end,1:gap:end,:);
    yD = yD-min(yD(:));
    zD = Z(1:gap:end,1:gap:end,:);
    zD = zD-min(zD(:));
    vxD = vx(1:gap:end,1:gap:end,:);
    vyD = vy(1:gap:end,1:gap:end,:);
    vzD = vz(1:gap:end,1:gap:end,:);

    %% Eigenvectors to Identify Axes

    vecs = regionprops3(relMask,'EigenVectors','EigenValues');
    if size(vecs,1) > 1
        error('too many objects')
    end
    vals = table2array(vecs(:,2));
    vals = vals{1};
    vecs = table2array(vecs(:,1));
    vecs = vecs{1};

    %% 3D Quiver - AP
    
    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));
    vD = vD(1:gap:end,1:gap:end,:); 

    figure(1)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD>=bins(jj)) & (vD<bins(jj+1));
        j = quiver3(xD(slice),yD(slice),zD(slice),vxD(slice),vyD(slice),vzD(slice),0,...
            'Color',cmap(jj,:),'LineWidth',0.25);
        if qscale ~=0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            hW = get(j,'WData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV,'WData',qscale*hW)
        end
        hold on
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'zdir','reverse')
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    set(gca,'GridColor',0.6*[1 1 1],'GridAlpha',0.5,'GridLineWidth',1)

    xlim([0 x2-x1]*xyscale)
    ylim([0 y2-y1]*xyscale)
    zlim([0 z2-z1]*zscale)
    set(gca,'XTick',0:50:1000,'YTick',0:50:1000,'ZTick',0:50:1000,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-80 80]) % Anterior = Left, Posterior = Right, Dorsal = Up, and Ventral = Down

    set(gcf,'Color','w')
    drawnow;

    quivAP = getframe(gca);
    camPos = get(gca,'CameraPosition');
    camTar = get(gca,'CameraTarget');
    camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end


    %% 3D Quiver - DV
    
    bins = linspace(-dvMax,dvMax,15);
    cmap = colorcet('D10'); % blue towards ventral, pink towards dorsal
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,2)) + (vx.*vecs(2,2)) + (vz.*vecs(3,2));
    vD = vD(1:gap:end,1:gap:end,:); 

    figure(2)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD>=bins(jj)) & (vD<bins(jj+1));
        j = quiver3(xD(slice),yD(slice),zD(slice),vxD(slice),vyD(slice),vzD(slice),0,...
            'Color',cmap(jj,:),'LineWidth',0.25);
        if qscale ~=0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            hW = get(j,'WData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV,'WData',qscale*hW)
        end
        hold on
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'zdir','reverse')
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    set(gca,'GridColor',0.6*[1 1 1],'GridAlpha',0.5,'GridLineWidth',1)

    xlim([0 x2-x1]*xyscale)
    ylim([0 y2-y1]*xyscale)
    zlim([0 z2-z1]*zscale)
    set(gca,'XTick',0:50:1000,'YTick',0:50:1000,'ZTick',0:50:1000,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-80 80]) % Anterior = Left, Posterior = Right, Dorsal = Up, and Ventral = Down

    set(gcf,'Color','w')
    drawnow;

    quivDV = getframe(gca);
    % camPos = get(gca,'CameraPosition');
    % camTar = get(gca,'CameraTarget');
    % camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end

    %% Crop Images for Display    

    im1plot = im1_ch0(y1:y2,x1:x2,z1:z2);  
    im1plot = imadjustn(im1plot,[100 600]/(2^16-1));
    im1plot = imgaussfilt(im1plot,0.05);

    imGrid = zeros(size(im1plot));
    cGrid = 20000;
    dxGrid = round(50/xyscale);
    dzGrid = round(50/zscale);

    imGrid(1:dxGrid:end,:,end) = cGrid;
    imGrid(:,1:dxGrid:end,end) = cGrid;
    imGrid(2:dxGrid:end,:,end) = cGrid;
    imGrid(:,2:dxGrid:end,end) = cGrid;
    imGrid(3:dxGrid:end,:,end) = cGrid;
    imGrid(:,3:dxGrid:end,end) = cGrid;

    imGrid(1:dxGrid:end,end,:) = cGrid;
    imGrid(:,end,1:dzGrid:end) = cGrid;
    imGrid(2:dxGrid:end,end,:) = cGrid;
    imGrid(:,end,2:dzGrid:end) = cGrid;
    imGrid(3:dxGrid:end,end,:) = cGrid;
    imGrid(:,end,3:dzGrid:end) = cGrid;

    imGrid(1,1:dxGrid:end,:) = cGrid;
    imGrid(1,:,1:dzGrid:end) = cGrid;
    imGrid(1,2:dxGrid:end,:) = cGrid;
    imGrid(1,:,2:dzGrid:end) = cGrid;
    imGrid(1,3:dxGrid:end,:) = cGrid;
    imGrid(1,:,3:dzGrid:end) = cGrid;

    imGrid(1,end,:) = cGrid;
    imGrid(1,:,end) = cGrid;
    imGrid(:,end,end) = cGrid;
    imGrid(end,:,end) = cGrid;
    imGrid(end,end,:) = cGrid;

    %%% Make 3D Rendering
    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [0 0 0];
    viewer.GradientColor = [0.2 0.2 0.2];
    viewer.Parent.Position = [100 100 size(quivAP.cdata,2) size(quivAP.cdata,1)];
    % viewer.Parent.Position = [100 100 1250 1020];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';

    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;

    viewer.CameraZoom = 1.78; % Emprically determined, not sure how to better address

    % pause(1) % this is to replicate "drawnow" but for a viewer
    volshow(im1plot,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling);

    pause(1)
    volshow(imGrid,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling,DisplayRange=[0 2^16]);  

    pause(2)
    frame = getframe(viewer.Parent);

    delete(viewer.Parent); % close the window now that the frame has been saved.

    %% Format 3D rendering window

    % Make an image of just the grid for cropping
    viewer = viewer3d;
    viewer.BackgroundColor = [0 0 0];
    viewer.BackgroundGradient = "off";
    viewer.Parent.Position = [100 100 size(quivAP.cdata,2) size(quivAP.cdata,1)];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';
    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;
    viewer.CameraZoom = 1.78; % Emprically determined, not sure how to better address
    volshow(imGrid,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling,DisplayRange=[0 2^16]);

    pause(2)
    frame2 = getframe(viewer.Parent);
    delete(viewer.Parent); % close the window now that the frame has been saved.

    frame2 = sum(frame2.cdata,3);
    frame2 = frame2>10;
    frame2 = imdilate(frame2,strel('disk',2));
    frame2 = imfill(frame2,'holes');

    frame2 = repmat(frame2,[1 1 3]);

    imPlotWhite = frame.cdata;
    imPlotWhite(~frame2) = 255;

    % figure(100)
    % imshow(imPlotWhite)   

    %% Make the movie panel
    
    strip = uint8(255*ones(size(imPlotWhite,1),50,3));
    panel = [imPlotWhite strip quivAP.cdata];

    % Add timestamp for movie
    tText = (exampleFrame-1)*tscale;
    mm = floor(tText);
    ss = tText-mm;
    ss = round(ss*60); % seconds
    tText = [num2str(mm,'%02u') ':' num2str(ss,'%02u')];

    a = size(panel);
    panel = insertText(panel,[a(2)-220 0],tText,BoxOpacity=0,TextColor='black',FontSize=60);

    % figure(3)
    % imshow(panel)

    imwrite(panel,[saveDir8 filesep 'MovieS8_AnteriorPosterior_frame' num2str(exampleFrame,'%03u') '.png'])

    %%
    strip = uint8(255*ones(size(imPlotWhite,1),50,3));
    panel = [imPlotWhite strip quivDV.cdata];

    % Add timestamp for movie
    tText = (exampleFrame-1)*tscale;
    mm = floor(tText);
    ss = tText-mm;
    ss = round(ss*60); % seconds
    tText = [num2str(mm,'%02u') ':' num2str(ss,'%02u')];

    a = size(panel);
    panel = insertText(panel,[a(2)-220 0],tText,BoxOpacity=0,TextColor='black',FontSize=60);

    % figure(3)
    % imshow(panel)

    imwrite(panel,[saveDir7 filesep 'MovieS7_DorsalVentral_frame' num2str(exampleFrame,'%03u') '.png'])


    %% Clean up

    close all
    disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end
