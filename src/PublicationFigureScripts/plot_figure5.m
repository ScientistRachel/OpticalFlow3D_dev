clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\SiMView_Drosophila';
dataName = 'fused_tp__ch_0';
dataDir = [imDir filesep 'OpticalFlow3D' filesep dataName];
saveDir = 'X:\Force Project\PublicationData\SiMView_Drosophila\PublicationFigures\MovieS7';
figDir = 'X:\Force Project\PublicationData\SiMView_Drosophila\PublicationFigures\Figure5';

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

if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
if ~exist(figDir,'dir')
    mkdir(figDir)
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
for exampleFrame = 54%4:88

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Load Images

    im1_ch0 = TIFFvolume([imDir filesep 'fused_tp_' num2str(exampleFrame-1) '_ch_0.tif'],Nz);
    im2_ch0 = TIFFvolume([imDir filesep 'fused_tp_' num2str(exampleFrame) '_ch_0.tif'],Nz);
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
    
    %% Useful quantities
    
    % mag = sqrt(vx.^2 + vy.^2 + vz.^2);   
    % theta = atan2(vy,vx);
    % phi = atan(vz./sqrt(vx.^2+vy.^2));

    % shared quiver values
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
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
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
%%
    quivAP = getframe(gca);
    imwrite(quivAP.cdata,[figDir filesep 'Figure5B_APflow_frame' num2str(exampleFrame,'%03u') '.png'])
    camPos = get(gca,'CameraPosition');
    camTar = get(gca,'CameraTarget');
    camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end

    %% 3D Quiver - DV
    
    bins = linspace(-dvMax,dvMax,15);
    cmap = colorcet('D10'); % blue towards ventral, pink towards dorsal
    % cDor = cmap(end,:);
    % cVen = cmap(1,:);
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
    imwrite(quivDV.cdata,[figDir filesep 'Figure5C_DVflow_frame' num2str(exampleFrame,'%03u') '.png'])
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

    imwrite(imPlotWhite,[figDir filesep 'Figure5A_Image_frame' num2str(exampleFrame,'%03u') '.png'])

    %% 3D Quiver - AP - Half - Z
    
    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));    

    xD1 = X(1:gap:end,1:gap:end,73:156);
    xD1 = xD1-min(xD1(:)); % keep the same zeroing as above to have same axes
    yD1 = Y(1:gap:end,1:gap:end,73:156);
    yD1 = yD1-min(yD1(:));
    zD1 = Z(1:gap:end,1:gap:end,73:156);
    zD1 = zD1-min(zD(:)); 
    vxD1 = vx(1:gap:end,1:gap:end,73:156);
    vyD1 = vy(1:gap:end,1:gap:end,73:156);
    vzD1 = vz(1:gap:end,1:gap:end,73:156);

    vD1 = vD(1:gap:end,1:gap:end,73:156); 

    figure(4)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
        j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
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
    % set(gca,'View',[-83 51]) % This view highlights the inside of the embryo

    set(gcf,'Color','w')
    drawnow;
%%
    quivAPhalf = getframe(gca);
    imwrite(quivAPhalf.cdata,[figDir filesep 'Figure5D_APflow_halfZ_frame' num2str(exampleFrame,'%03u') '.png'])
    % camPos = get(gca,'CameraPosition');
    % camTar = get(gca,'CameraTarget');
    % camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end

 %% 3D Quiver - AP - Slice - Z
    
    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));   

    % z1Slice = 70;
    % z2Slice = 110;

    z1Slice = 67;
    z2Slice = 107;

    xD1 = X(1:gap:end,1:gap:end,z1Slice:z2Slice);
    xD1 = xD1-min(xD1(:)); % keep the same zeroing as above to have same axes
    yD1 = Y(1:gap:end,1:gap:end,z1Slice:z2Slice);
    yD1 = yD1-min(yD1(:));
    zD1 = Z(1:gap:end,1:gap:end,z1Slice:z2Slice);
    zD1 = zD1-min(zD(:)); 
    vxD1 = vx(1:gap:end,1:gap:end,z1Slice:z2Slice);
    vyD1 = vy(1:gap:end,1:gap:end,z1Slice:z2Slice);
    vzD1 = vz(1:gap:end,1:gap:end,z1Slice:z2Slice);

    vD1 = vD(1:gap:end,1:gap:end,z1Slice:z2Slice); 

    figure(4)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
        j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
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
    % set(gca,'View',[-83 51]) % This view highlights the inside of the embryo

    set(gcf,'Color','k')
    drawnow;

    quivAPslab = getframe(gca);
    imwrite(quivAPslab.cdata,[figDir filesep 'Figure5D_APflow_Zslab_frame' num2str(exampleFrame,'%03u') '.png'])
    % camPos = get(gca,'CameraPosition');
    % camTar = get(gca,'CameraTarget');
    % camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end


    %% 3D Quiver - DV - Half
% 
%     bins = linspace(-dvMax,dvMax,15);
%     cmap = colorcet('D10'); % blue towards ventral, pink towards dorsal
%     % cDor = cmap(end,:);
%     % cVen = cmap(1,:);
%     cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);
% 
%     vD = (vy.*vecs(1,2)) + (vx.*vecs(2,2)) + (vz.*vecs(3,2));    
% 
%     xD1 = X(1:gap:end,1:gap:end,73:156);
%     xD1 = xD1-min(xD1(:)); % keep the same zeroing as above to have same axes
%     yD1 = Y(1:gap:end,1:gap:end,73:156);
%     yD1 = yD1-min(yD1(:));
%     zD1 = Z(1:gap:end,1:gap:end,73:156);
%     zD1 = zD1-min(zD(:)); 
%     vxD1 = vx(1:gap:end,1:gap:end,73:156);
%     vyD1 = vy(1:gap:end,1:gap:end,73:156);
%     vzD1 = vz(1:gap:end,1:gap:end,73:156);
% 
%     vD1 = vD(1:gap:end,1:gap:end,73:156); 
% 
%     figure(5)
%     set(gcf,'Position',[300 150 1250 1020])
%     for jj = 1:length(bins)-1     
%         slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
%         j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
%             'Color',cmap(jj,:),'LineWidth',0.25);
%         if qscale ~=0
%             hU = get(j,'UData') ;
%             hV = get(j,'VData') ;
%             hW = get(j,'WData') ;
%             set(j,'UData',qscale*hU,'VData',qscale*hV,'WData',qscale*hW)
%         end
%         hold on
%     end
%     hold off
%     set(gca,'DataAspectRatio',[1 1 1])
%     set(gca,'zdir','reverse')
%     set(gca,'ydir','reverse')
%     set(gca,'Color','k')
%     set(gca,'GridColor',0.6*[1 1 1],'GridAlpha',0.5,'GridLineWidth',1)
% 
%     xlim([0 x2-x1]*xyscale)
%     ylim([0 y2-y1]*xyscale)
%     zlim([0 z2-z1]*zscale)
%     set(gca,'XTick',0:50:1000,'YTick',0:50:1000,'ZTick',0:50:1000,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
%     set(gca,'View',[-80 80]) % Anterior = Left, Posterior = Right, Dorsal = Up, and Ventral = Down
%     % set(gca,'View',[-83 51]) % This view highlights the inside of the embryo
% 
%     set(gcf,'Color','w')
%     drawnow;
% %%
%     quivDVhalf = getframe(gca);
%     imwrite(quivDVhalf.cdata,[figDir filesep 'Figure5D_DVflow_halfZ_frame' num2str(exampleFrame,'%03u') '.png'])
%     % camPos = get(gca,'CameraPosition');
%     % camTar = get(gca,'CameraTarget');
%     % camUp = get(gca,'CameraUpVector');
% 
%     if cleanFig
%         close all
%     end

%% 3D Quiver - AP - Half - Y

    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));    

    ycrop = 750;

    xD1 = X(1:gap:ycrop,1:gap:end,:);
    xD1 = xD1-min(xD1(:)); % keep the same zeroing as above to have same axes
    yD1 = Y(1:gap:ycrop,1:gap:end,:);
    yD1 = yD1-min(Y(:));
    zD1 = Z(1:gap:ycrop,1:gap:end,:);
    zD1 = zD1-min(zD1(:)); 
    vxD1 = vx(1:gap:ycrop,1:gap:end,:);
    vyD1 = vy(1:gap:ycrop,1:gap:end,:);
    vzD1 = vz(1:gap:ycrop,1:gap:end,:);

    vD1 = vD(1:gap:ycrop,1:gap:end,:); 

    figure(500)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
        j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
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
    % set(gca,'View',[-83 51]) % This view highlights the inside of the embryo

    set(gcf,'Color','w')
    drawnow;
% %%
%     quivAPhalf = getframe(gca);
%     imwrite(quivAPhalf.cdata,[figDir filesep 'Figure5D_APflow_halfZ_frame' num2str(exampleFrame,'%03u') '.png'])
%     % camPos = get(gca,'CameraPosition');
%     % camTar = get(gca,'CameraTarget');
%     % camUp = get(gca,'CameraUpVector');
% 
%     if cleanFig
%         close all
%     end

%% 3D Quiver - AP - Half - Dorsal

    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));    

    xcrop = 225;

    xD1 = X(1:gap:end,xcrop:gap:end,:);
    xD1 = xD1-min(X(:)); % keep the same zeroing as above to have same axes
    yD1 = Y(1:gap:end,xcrop:gap:end,:);
    yD1 = yD1-min(Y(:));
    zD1 = Z(1:gap:end,xcrop:gap:end,:);
    zD1 = zD1-min(Z(:)); 
    vxD1 = vx(1:gap:end,xcrop:gap:end,:);
    vyD1 = vy(1:gap:end,xcrop:gap:end,:);
    vzD1 = vz(1:gap:end,xcrop:gap:end,:);

    vD1 = vD(1:gap:end,xcrop:gap:end,:); 

    figure(500)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
        j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
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
    % set(gca,'View',[-83 51]) % This view highlights the inside of the embryo

    set(gcf,'Color','w')
    drawnow;
% %%
%     quivAPhalf = getframe(gca);
%     imwrite(quivAPhalf.cdata,[figDir filesep 'Figure5D_APflow_halfZ_frame' num2str(exampleFrame,'%03u') '.png'])
%     % camPos = get(gca,'CameraPosition');
%     % camTar = get(gca,'CameraTarget');
%     % camUp = get(gca,'CameraUpVector');
% 
%     if cleanFig
%         close all
%     end

    %% Make 3D Rendering of relmask - Full Fly

    cartoon = imgaussfilt(double(relMask(:,:,:)),6);

    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [1 1 1];
    viewer.GradientColor = [1 1 1];
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
    volshow(cartoon,Parent=viewer,Transformation=scaling,Colormap=bone(255));

    pause(2)
    frame = getframe(viewer.Parent);

    delete(viewer.Parent); % close the window now that the frame has been saved.

    imwrite(frame.cdata,[figDir filesep 'Figure5_legend_fullEmbryo_frame' num2str(exampleFrame,'%03u') '.png'])

    %% Half Fly

    cartoon = imgaussfilt(double(relMask(:,:,73:156)),6);

    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [1 1 1];
    viewer.GradientColor = [1 1 1];
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
    volshow(cartoon,Parent=viewer,Transformation=scaling,Colormap=bone(255));

    pause(2)
    frame = getframe(viewer.Parent);

    delete(viewer.Parent); % close the window now that the frame has been saved.

    imwrite(frame.cdata,[figDir filesep 'Figure5_legend_halfEmbryo_frame' num2str(exampleFrame,'%03u') '.png'])

    %% Slab Fly

    cartoon = imgaussfilt(double(relMask(:,:,z1Slice:z2Slice)),6);

    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [1 1 1];
    viewer.GradientColor = [1 1 1];
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
    volshow(cartoon,Parent=viewer,Transformation=scaling,Colormap=bone(255));

    pause(2)
    frame = getframe(viewer.Parent);

    delete(viewer.Parent); % close the window now that the frame has been saved.

    imwrite(frame.cdata,[figDir filesep 'Figure5_legend_halfEmbryo_frame' num2str(exampleFrame,'%03u') '.png'])

    %% Outline the fly on the half figure
    % 
    % bins = linspace(-apMax,apMax,15);
    % cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);
    % 
    % vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));
    % vD = vD(1:gap:end,1:gap:end,:); 
    % 
    % figure(1)
    % set(gcf,'Position',[300 150 1250 1020])
    % for jj = 1:length(bins)-1     
    %     slice = (vD>=bins(jj)) & (vD<bins(jj+1));
    %     j = quiver3(xD(slice),yD(slice),zD(slice),vxD(slice),vyD(slice),vzD(slice),0,...
    %         'Color',cmap(jj,:),'LineWidth',0.25);
    %     if qscale ~=0
    %         hU = get(j,'UData') ;
    %         hV = get(j,'VData') ;
    %         hW = get(j,'WData') ;
    %         set(j,'UData',qscale*hU,'VData',qscale*hV,'WData',qscale*hW)
    %     end
    %     hold on
    % end
    % hold off
    % set(gca,'DataAspectRatio',[1 1 1])
    % set(gca,'zdir','reverse')
    % set(gca,'ydir','reverse')
    % set(gca,'Color','k')
    % set(gca,'GridColor',0.6*[1 1 1],'GridAlpha',0.5,'GridLineWidth',1)
    % 
    % xlim([0 x2-x1]*xyscale)
    % ylim([0 y2-y1]*xyscale)
    % zlim([0 z2-z1]*zscale)
    % axis off
    % set(gca,'View',[-80 80]) % Anterior = Left, Posterior = Right, Dorsal = Up, and Ventral = Down
    % 
    % set(gcf,'Color','k')
    % drawnow;
    % 
    % outIm = getframe(gca);
    % close all
    % 
    % outline = outIm.cdata;
    % outline = outline>0;
    % outline = sum(outline,3);
    % B = bwboundaries(outline);
    % imshow(quivAPhalf.cdata)
    % hold on
    % plot(B{1}(:,2),B{1}(:,1),'--','Color',[240, 228, 66]/255,'LineWidth',3)
    % hold off




    %% Inset
    % v0: near the back wall
    % ylim([290 350])
    % xlim([200 201])
    % view([90 0])
    % zlim([110 150])

    % v1: nice swirl, not sure about hte image
    % % xI1 = 320;
    % % xI2 = 369;
    % yI1 = 1057;
    % yI2 = 1106;
    % zI1 = 127;
    % zI2 = 131;

    xI1 = 345;
    xI2 = 375;
    yI1 = 1030;
    yI2 = 1060;
    zI1 = 120;
    zI2 = 130;
    qscale1 = 0.5;

    % v1: near where the inset figures would be anyway

    bins = linspace(-apMax,apMax,15);
    cmap = colorcet('D13'); % Green points towards anterior, blue towards posterior
    % cAnt = cmap(end,:);
    % cPost = cmap(1,:);
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = (vy.*vecs(1,1)) + (vx.*vecs(2,1)) + (vz.*vecs(3,1));    

    xD1 = X(yI1:yI2,xI1:xI2,zI1:zI2);
    xD1 = xD1-min(xD1(:));
    yD1 = Y(yI1:yI2,xI1:xI2,zI1:zI2);
    yD1 = yD1-min(yD1(:));
    zD1 = Z(yI1:yI2,xI1:xI2,zI1:zI2);
    zD1 = zD1-min(zD(:)); 
    vxD1 = vx(yI1:yI2,xI1:xI2,zI1:zI2);
    vyD1 = vy(yI1:yI2,xI1:xI2,zI1:zI2);
    vzD1 = vz(yI1:yI2,xI1:xI2,zI1:zI2);

    vD1 = vD(yI1:yI2,xI1:xI2,zI1:zI2);

    figure(4)
    set(gcf,'Position',[300 150 1250 1020])
    for jj = 1:length(bins)-1     
        slice = (vD1>=bins(jj)) & (vD1<bins(jj+1));
        j = quiver3(xD1(slice),yD1(slice),zD1(slice),vxD1(slice),vyD1(slice),vzD1(slice),0,...
            'Color',cmap(jj,:),'LineWidth',1.5);
        if qscale1 ~=0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            hW = get(j,'WData') ;
            set(j,'UData',qscale1*hU,'VData',qscale1*hV,'WData',qscale1*hW)
        end
        hold on
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'zdir','reverse')
    set(gca,'xdir','reverse')
    set(gca,'ydir','normal')
    set(gca,'Color','k')
    set(gca,'GridColor',0.6*[1 1 1],'GridAlpha',0.5,'GridLineWidth',1)
    % 
    % ylim([-1 22])
    % xlim([-1 21])
    % set(gca,'XTick',0:10:100,'YTick',0:10:40,'ZTick',0:10:200,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[90 90])

    set(gcf,'Color','w')
    drawnow;
% %%
%     quivAPinset = getframe(gca);
%     imwrite(quivAPint.cdata,[figDir filesep 'Figure5D_APflow_inset_frame' num2str(exampleFrame,'%03u') '.png'])
%     % camPos = get(gca,'CameraPosition');
%     % camTar = get(gca,'CameraTarget');
%     % camUp = get(gca,'CameraUpVector');
% 
%     if cleanFig
%         close all
%     end
%

    testIm = im1_ch0(yI1:yI2,xI1:xI2,zI1:zI2);
    testIm2 = im2_ch0(yI1:yI2,xI1:xI2,zI1:zI2);

    testIm = max(testIm,[],3);
    testIm = imresize(testIm,4,'nearest');
    testIm2 = max(testIm2,[],3);
    testIm2 = imresize(testIm2,4,'nearest');

    figure(5)
    imshowpair(imadjust(testIm),imadjust(testIm2))
    set(gca,'View',[90 90])
    set(gca,'xdir','reverse')
    set(gca,'ydir','normal')

    

    %% Clean up

    close all
    disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end

%% Make some colorbars
% 
% figure(7)
% set(gcf,'Position',[695   658   754   420])
% imagesc(rand(10))
% clim(apMax*[-1 1])
% colormap(colorcet('D13'))
% h = colorbar;
% set(get(h,'Label'),'String','AP-axis Magnitude')
% set(gca,'FontSize',20)
% saveas(gcf,[figDir filesep 'AP_scaleBar.png'])
% saveas(gcf,[figDir filesep 'AP_scaleBar.svg'])
% 
% figure(8)
% set(gcf,'Position',[695   658   500 500])
% imagesc(rand(10))
% clim(dvMax*[-1 1])
% colormap(colorcet('D10'))
% h = colorbar('Location','SouthOutside');
% set(get(h,'Label'),'String','DV-axis Magnitude')
% set(gca,'FontSize',20)
% saveas(gcf,[figDir filesep 'DV_scaleBar.png'])
% saveas(gcf,[figDir filesep 'DV_scaleBar.svg'])
% 
% close all
% disp('script complete')