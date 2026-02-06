clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\deskew_after_decon\fused';
dataName = 'RotatedCropped_fused_tp__ch_0';
dataDir = [imDir filesep 'OpticalFlow3D_MATLAB' filesep dataName];
saveDir5 = 'X:\Force Project\PublicationData\MOSAIC_Actin\PublicationFigures\MovieS5';
saveDir6 = 'X:\Force Project\PublicationData\MOSAIC_Actin\PublicationFigures\MovieS6';
figDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\PublicationFigures\Figure4';

% Metadata
xyscale = 0.108; % um/pixel
zscale = 0.5*sin(32.45*pi/180);
tscale = 60075/1000/60; % minutes/frame

% reliability thresholding
relPer = 95;
minSize = 10^5;

% For plots
qscale = 1; % vector scaling, 1 or 0 both equal actual displacement
gap = 1; % show every vector or only every gap'th vector?
magMax = 1.25; % um/min

% Large ROI version
x1 = 1050;
x2 = 1600;
y1 = 550;
y2 = 1000;

x3 = 325;
x4 = 1050;
y3 = 500;
y4 = 1250;

z1 = 15;
z2 = 113;

keyFrames = [4, 25, 45, 73]; % These are the frames called out in figures, not just the movies

% clean figures
cleanFig = true;

%% Set up the save directories

if ~exist(saveDir5,'dir')
    mkdir(saveDir5)
end
if ~exist(saveDir6,'dir')
    mkdir(saveDir6)
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

%% Loop through frames
for exampleFrame = 4:75

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Load Images

    im1 = TIFFvolume([imDir filesep 'RotatedCropped_fused_tp_' num2str(exampleFrame-1,'%02u') '_ch_0.tif'],Nz);
    im1N = TIFFvolume([imDir filesep 'RotatedCropped_fused_tp_' num2str(exampleFrame-1,'%02u') '_ch_1.tif'],Nz);
    % disp('raw images loaded')

    %% Reliability Thresholding
    
    relFile = [dataName '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    rel = TIFFvolume([dataDir filesep relFile],Nz);
    % disp('rel loaded')

    relThresh = prctile(rel(:),relPer);
    % disp(['Reliability Threshold: ' num2str(relThresh)])
    
    relMask = rel > relThresh;
    relMask = bwareaopen(relMask,minSize);
    
    %% Load in velocities
    
    vx = TIFFvolume([dataDir filesep dataName '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx = vx.*relMask./relMask*xyscale/tscale;
    % disp('vx loaded')
    
    vy = TIFFvolume([dataDir filesep dataName '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy = vy.*relMask./relMask*xyscale/tscale;
    % disp('vy loaded')
    
    vz = TIFFvolume([dataDir filesep dataName '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vz = vz.*relMask./relMask*zscale/tscale;
    % disp('vz loaded')
    
    %% Useful quantities
    
    mag = sqrt(vx.^2 + vy.^2 + vz.^2);   
    theta = atan2(vy,vx);
    phi = atan(vz./sqrt(vx.^2+vy.^2));

    %% Dividing Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% 3D Quiver - Dividing Cell

    % shared quiver values
    xD = X(y1:gap:y2,x1:gap:x2,z1:gap:z2);
    xD = xD-min(xD(:)); % Zeroing here lets me match the quiver and image figure camera positions effectively.
    yD = Y(y1:gap:y2,x1:gap:x2,z1:gap:z2);
    yD = yD-min(yD(:));
    zD = Z(y1:gap:y2,x1:gap:x2,z1:gap:z2);
    zD = zD-min(zD(:));
    vxD = vx(y1:gap:y2,x1:gap:x2,z1:gap:z2);
    vyD = vy(y1:gap:y2,x1:gap:x2,z1:gap:z2);
    vzD = vz(y1:gap:y2,x1:gap:x2,z1:gap:z2);

    % 3D Quiver - Dividing Cell - Magnitude
    bins = linspace(0,magMax,15-1);
    bins(15) = Inf;
    cmap = colorcet('L8');
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = mag(y1:gap:y2,x1:gap:x2,z1:gap:z2);    

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

    xlim([0 60]) % Round to 60 for the x limits so that the grid lines make sense (from 59.4 --> 60)
    ylim([0 y2-y1]*xyscale)
    zlim([0 z2-z1]*zscale)
    set(gca,'XTick',0:5:65,'YTick',0:5:55,'ZTick',0:5:25,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-37 53.5])

    set(gcf,'Color','w')
    drawnow;

    quivMag = getframe(gca);
    camPos = get(gca,'CameraPosition');
    camTar = get(gca,'CameraTarget');
    camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end

    %% 3D Quiver - Dividing Cell - Phi
    
    bins = linspace(-pi/2-eps,pi/2+eps,15);
    cmap = colorcet('D2');
    cmap = flipud(cmap); % magenta up
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = phi(y1:y2,x1:x2,z1:z2);

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
    
    xlim([0 60]) % Round to 60 for the x limits so that the grid lines make sense (from 59.4 --> 60)
    ylim([0 y2-y1]*xyscale)
    zlim([0 z2-z1]*zscale)
    set(gca,'XTick',0:5:65,'YTick',0:5:55,'ZTick',0:5:25,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-37 53.5])
    
    set(gcf,'Color','w')
    drawnow;

    quivPhi = getframe(gca);

    if cleanFig
        close all
    end

    %% Crop Images for Display    

    imNplot = im1N(y1:y2,x1:x2,z1:z2);
    imNplot = imadjustn(imNplot,[15 75]/(2^16-1));
    imNplot = imgaussfilt(imNplot,1);

    im1plot = im1(y1:y2,x1:x2,z1:z2);  
    im1plot = imadjustn(im1plot,[5 420]/(2^16-1));
    im1plot = imgaussfilt(im1plot,0.05);

    imGrid = zeros(size(im1plot));
    cGrid = 20000;
    dxGrid = round(5/xyscale);
    dzGrid = round(5/zscale);

    imGrid(1:dxGrid:end,:,end) = cGrid;
    imGrid(:,1:dxGrid:end,end) = cGrid;

    imGrid(1:dxGrid:end,end,5:end) = cGrid;
    imGrid(:,end,end:-dzGrid:1) = cGrid;

    imGrid(1,1:dxGrid:end,5:end) = cGrid;
    imGrid(1,:,end:-dzGrid:1) = cGrid;

    imGrid(1,end,5:end) = cGrid;
    imGrid(1,:,end) = cGrid;
    imGrid(:,end,end) = cGrid;
    imGrid(end,:,end) = cGrid;
    imGrid(end,end,5:end) = cGrid;

    %% Make 3D Rendering
    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [0 0 0];
    viewer.GradientColor = [0.2 0.2 0.2];
    viewer.Parent.Position = [100 100 size(quivMag.cdata,2) size(quivMag.cdata,1)];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';

    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;

    viewer.CameraZoom = 1.07; % Emprically determined, not sure how to better address

    % pause(1)
    volshow(imNplot,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Colormap=abyss(256),Transformation=scaling);

    pause(1) % this is to replicate "drawnow" but for a viewer
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
    viewer.Parent.Position = [100 100 size(quivMag.cdata,2) size(quivMag.cdata,1)];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';
    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;
    viewer.CameraZoom = 1.07; % Emprically determined, not sure how to better address
    volshow(imGrid,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling,DisplayRange=[0 2^16]);

    pause(2)
    frame2 = getframe(viewer.Parent);
    delete(viewer.Parent); % close the window now that the frame has been saved.

    frame2 = sum(frame2.cdata,3);
    frame2 = frame2>200;
    frame2 = imdilate(frame2,strel('disk',2));
    frame2 = imfill(frame2,'holes');
    frame2 = repmat(frame2,[1 1 3]);

    imPlotWhite = frame.cdata;
    imPlotWhite(~frame2) = 255;

    % figure;
    % imshow(imPlotWhite)

    if sum(exampleFrame==keyFrames)
        imwrite(imPlotWhite,[figDir filesep 'Fig3A_frame' num2str(exampleFrame,'%03u') '.png'])
    end

    %% Make the movie panel
    
    strip = uint8(255*ones(size(imPlotWhite,1),50,3));
    panel = [imPlotWhite strip quivMag.cdata strip quivPhi.cdata];

    % Add timestamp for movie
    tText = (exampleFrame-1)*tscale;
    hh = floor(tText/60);
    mm = tText - 60*hh;
    mm = floor(mm);

    tText = [num2str(hh,'%02u') ':' num2str(mm,'%02u')];
    panel = insertText(panel,[20 15],tText,BoxOpacity=0,TextColor='black',FontSize=60);
 
    % figure(3)
    % imshow(panel)

    imwrite(panel,[saveDir5 filesep 'MovieS5_frame' num2str(exampleFrame,'%03u') '.png'])


    %% Migrating Cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% 3D Quiver - Migrating Cell

    % shared quiver values
    xD = X(y3:gap:y4,x3:gap:x4,z1:gap:z2);
    xD = xD-min(xD(:)); % Zeroing here lets me match the quiver and image figure camera positions effectively.
    yD = Y(y3:gap:y4,x3:gap:x4,z1:gap:z2);
    yD = yD-min(yD(:));
    zD = Z(y3:gap:y4,x3:gap:x4,z1:gap:z2);
    zD = zD-min(zD(:));
    vxD = vx(y3:gap:y4,x3:gap:x4,z1:gap:z2);
    vyD = vy(y3:gap:y4,x3:gap:x4,z1:gap:z2);
    vzD = vz(y3:gap:y4,x3:gap:x4,z1:gap:z2);

    % 3D Quiver - Magnitude
    bins = linspace(0,magMax,15-1);
    bins(15) = Inf;
    cmap = colorcet('L8');
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = mag(y3:gap:y4,x3:gap:x4,z1:gap:z2);    

    figure(4)
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
    
    xlim([0 x4-x3]*xyscale)
    ylim([0 y4-y3]*xyscale)
    zlim([0 z2-z1]*zscale)
    ztick = ((z2-z1)*zscale):-5:0;
    ztick = sort(ztick);
    set(gca,'XTick',0:5:100,'YTick',0:5:100,'ZTick',ztick,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-35 50])

    set(gcf,'Color','w')
    drawnow;

    quivMag = getframe(gca);
    camPos = get(gca,'CameraPosition');
    camTar = get(gca,'CameraTarget');
    camUp = get(gca,'CameraUpVector');

    if cleanFig
        close all
    end

    %% 3D Quiver - Migrating - Theta
    
    bins = linspace(-pi-eps,pi+eps,15);
    cmap = colorcet('C8');
    cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

    vD = theta(y3:y4,x3:x4,z1:z2);

    figure(5)
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
    
    xlim([0 x4-x3]*xyscale)
    ylim([0 y4-y3]*xyscale)
    zlim([0 z2-z1]*zscale)
    ztick = ((z2-z1)*zscale):-5:0;
    ztick = sort(ztick);
    set(gca,'XTick',0:5:100,'YTick',0:5:100,'ZTick',ztick,'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{},'TickLength',[0 0])
    set(gca,'View',[-35 50])
    
    set(gcf,'Color','w')
    drawnow;

    quivPhi = getframe(gca);

    if cleanFig
        close all
    end

    %% Crop Images for Display    

    imNplot = im1N(y3:y4,x3:x4,z1:z2);
    imNplot = imadjustn(imNplot,[15 75]/(2^16-1));
    imNplot = imgaussfilt(imNplot,1);

    im1plot = im1(y3:y4,x3:x4,z1:z2);  
    % im1plot = imadjustn(im1plot,[5 420]/(2^16-1));
    im1plot = imadjustn(im1plot,[5 250]/(2^16-1));
    im1plot = imgaussfilt(im1plot,0.05);

    imGrid = zeros(size(im1plot));
    cGrid = 20000;
    dxGrid = round(5/xyscale);
    dzGrid = round(5/zscale);

    imGrid(1:dxGrid:end,:,end) = cGrid;
    imGrid(:,1:dxGrid:end,end) = cGrid;

    imGrid(1:dxGrid:end,end,5:end) = cGrid;
    imGrid(:,end,end:-dzGrid:1) = cGrid;

    imGrid(1,1:dxGrid:end,5:end) = cGrid;
    imGrid(1,:,end:-dzGrid:1) = cGrid;

    imGrid(1,end,5:end) = cGrid;
    imGrid(1,:,end) = cGrid;
    imGrid(:,end,end) = cGrid;
    imGrid(end,:,end) = cGrid;
    imGrid(end,end,5:end) = cGrid;

    %% Make 3D Rendering
    viewer = viewer3d;

    scaling = affinetform3d([xyscale 0 0 0 ; 0 xyscale 0 0 ; 0 0 zscale 0 ; 0 0 0 1]);       

    viewer.BackgroundColor = [0 0 0];
    viewer.GradientColor = [0.2 0.2 0.2];
    viewer.Parent.Position = [200 200 size(quivMag.cdata,2) size(quivMag.cdata,1)];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';

    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;

    viewer.CameraZoom = 1.125; % Emprically determined, not sure how to better address

    % pause(1)
    volshow(imNplot,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Colormap=abyss(256),Transformation=scaling);

    pause(1) % this is to replicate "drawnow" but for a viewer
    volshow(im1plot,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling);

    pause(1)
    volshow(imGrid,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling,DisplayRange=[0 2^16]);  

    pause(1)

    frame = getframe(viewer.Parent);

    delete(viewer.Parent); % close the window now that the frame has been saved.

    %% Format 3D rendering window

    % Make an image of just the grid for cropping
    viewer = viewer3d;
    viewer.BackgroundColor = [0 0 0];
    viewer.BackgroundGradient = "off";
    viewer.Parent.Position = [200 200 size(quivMag.cdata,2) size(quivMag.cdata,1)];
    viewer.RenderingQuality = "high";
    viewer.OrientationAxes = 'off';
    viewer.CameraPosition = camPos;
    viewer.CameraPositionMode = 'auto';
    viewer.CameraTarget =  camTar;
    viewer.CameraTargetMode = 'auto';
    viewer.CameraUpVector = camUp;
    viewer.CameraZoom = 1.125; % Emprically determined, not sure how to better address
    volshow(imGrid,RenderingStyle="MaximumIntensityProjection",Parent=viewer,Transformation=scaling,DisplayRange=[0 2^16]);

    pause(1)
    frame2 = getframe(viewer.Parent);
    delete(viewer.Parent); % close the window now that the frame has been saved.

    frame2 = sum(frame2.cdata,3);
    frame2 = frame2>100;
    frame2 = imdilate(frame2,strel('disk',2));
    frame2 = imfill(frame2,'holes');
    frame2 = repmat(frame2,[1 1 3]);

    imPlotWhite = frame.cdata;
    imPlotWhite(~frame2) = 255;

    % figure(20)
    % imshow(imPlotWhite)

    if sum(exampleFrame==keyFrames)
        imwrite(imPlotWhite,[figDir filesep 'Fig3D_frame' num2str(exampleFrame,'%03u') '.png'])
    end

    %% Make the movie panel
    
    strip = uint8(255*ones(size(imPlotWhite,1),100,3));
    panel = [imPlotWhite strip quivMag.cdata strip quivPhi.cdata];

    panel = insertText(panel,[20 15],tText,BoxOpacity=0,TextColor='black',FontSize=60);

    % figure(6)
    % imshow(panel)

    imwrite(panel,[saveDir6 filesep 'MovieS6_frame' num2str(exampleFrame,'%03u') '.png'])
   
   
    %% Clean up

    close all
    disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end

%% Make some colorbars

figure(7)
set(gcf,'Position',[695   658   754   420])
imagesc(rand(10))
clim([0 magMax])
colormap(colorcet('L8'))
h = colorbar;
set(get(h,'Label'),'String',['Magnitude (' char(181) 'm/min)'])
set(gca,'FontSize',20)
saveas(gcf,[figDir filesep 'MagnitudeScaleBar.png'])
saveas(gcf,[figDir filesep 'MagnitudeScaleBar.svg'])

figure(8)
set(gcf,'Position',[695   658   754   420])
imagesc(rand(10))
clim([-pi pi]/2*180/pi)
colormap(colorcet('D2'))
h = colorbar;
h.Ticks = -90:45:90;
set(get(h,'Label'),'String','\phi (degrees)')
set(gca,'FontSize',20)
saveas(gcf,[figDir filesep 'PhiScaleBar.png'])
saveas(gcf,[figDir filesep 'PhiScaleBar.svg'])

close all
disp('script complete')