clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\deskew_after_decon\fused';
dataName = 'RotatedCropped_fused_tp__ch_0';
dataDir = [imDir filesep 'OpticalFlow3D_MATLAB' filesep dataName];
saveDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\PublicationFigures\Figure4';

% Metadata
xyscale = 0.108; % um/pixel
zscale = 0.5*sin(32.45*pi/180);
tscale = 60075/1000/60; % minutes/frame

% reliability thresholding
relPer = 95;
minSize = 10^5;

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

% Phi Binning
Nbins = 30;
phiBins = linspace(-pi/2,pi/2,Nbins+1); % Edges
phiCenters = phiBins(1:end-1) + phiBins(2)-phiBins(1); % Centers

% Preallocate storage
meanmag = NaN*ones(2,78);
distPhiW = NaN*ones(length(phiBins)-1,78,2);

if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

%% Set up general matrices rather than recalculating each time
relFile = [dataName '_rel_t' num2str(4,'%04u') '.tiff'];
meta = imfinfo([dataDir filesep relFile]);
Nx = meta(1).Width;
Ny = meta(2).Height;
Nz = length(meta);


%% Loop through frames
for exampleFrame = 4:75

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

    %% Reliability Thresholding
    
    relFile = [dataName '_rel_t' num2str(exampleFrame,'%04u') '.tiff'];
    rel = TIFFvolume([dataDir filesep relFile],Nz);
    % disp('rel loaded')

    relThresh = prctile(rel(:),relPer);

    % Crop after the threshold calc, for consistency with other figures
    % Cropping now saves memory for the calculations later    
    rel1 = rel(y3:y4,x3:x4,z1:z2); % migrating cell
    rel = rel(y1:y2,x1:x2,z1:z2); % dividing cell

    relMask = rel > relThresh;
    relMask = bwareaopen(relMask,minSize);

    relMask1 = rel1 > relThresh;
    relMask1 = bwareaopen(relMask1,minSize);
    
    
    %% Load in velocities

    % vx
    
    vx = TIFFvolume([dataDir filesep dataName '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);

    vx1 = vx(y3:y4,x3:x4,z1:z2);
    vx1 = vx1.*relMask1./relMask1*xyscale/tscale;

    vx = vx(y1:y2,x1:x2,z1:z2);
    vx = vx.*relMask./relMask*xyscale/tscale;

    % vy
    
    vy = TIFFvolume([dataDir filesep dataName '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);

    vy1 = vy(y3:y4,x3:x4,z1:z2);
    vy1 = vy1.*relMask1./relMask1*xyscale/tscale;

    vy = vy(y1:y2,x1:x2,z1:z2);
    vy = vy.*relMask./relMask*xyscale/tscale;

    % vz
    
    vz = TIFFvolume([dataDir filesep dataName '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);

    vz1 = vz(y3:y4,x3:x4,z1:z2);
    vz1 = vz1.*relMask1./relMask1*zscale/tscale;

    vz = vz(y1:y2,x1:x2,z1:z2);
    vz = vz.*relMask./relMask*zscale/tscale;

    %%% Directionality   
    phi = atan(vz./sqrt(vx.^2+vy.^2));  
    phi1 = atan(vz1./sqrt(vx1.^2+vy1.^2));
    
    %% Magnitude Over Time
    
    mag = sqrt(vx.^2 + vy.^2 + vz.^2);
    badRows = isnan(mag(:));
    meanmag(1,exampleFrame) = mean(mag(~badRows));

    mag1 = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
    badRows1 = isnan(mag1(:));
    meanmag(2,exampleFrame) = mean(mag1(~badRows1));

    %% Weighted phi

    subs = discretize(-phi(:),phiBins);
    histwPhi = accumarray(subs(~badRows),mag(~badRows),[length(phiBins),1]);    
    phiCounts = histwPhi(1:end-1)/length(mag(~badRows));    
    distPhiW(:,exampleFrame,1) = phiCounts;   

    subs = discretize(-phi1(:),phiBins);
    histwPhi = accumarray(subs(~badRows1),mag1(~badRows1),[length(phiBins),1]);    
    phiCounts = histwPhi(1:end-1)/length(mag1(~badRows1));    
    distPhiW(:,exampleFrame,2) = phiCounts;    
    
end

%% Save data for easier iteration on the plots

save([saveDir filesep 'figure4_timeTraces.mat'],'meanmag','distPhiW',...
    'xyscale','zscale','tscale',...
    'relPer','minSize',...
    'x1','x2','x3','x4','y1','y2','y3','y4',...
    'Nbins','phiBins','phiCenters')

disp([char(datetime('now')) ' - Data Saved'])


%% Mean Magnitude vs Time

t = (4:75)*tscale;

figure(1)
set(gcf,'Position',[400 500 560*1.5 420])
plot(t,meanmag(1,4:75),'Color',0.6*[1 1 1],'LineWidth',0.5)
hold on
scatter(t,meanmag(1,4:75),50,meanmag(1,4:75),'filled')
hold off
colormap(colorcet('L8'))
clim([0.15 0.4]) % This is different than the quiver graphs! The mean will be lower because there are always regions that are not moving.
ylabel(['Mean Magnitude (' char(181) 'm/min)'])
xlabel('Time (min)')
set(gca,'FontSize',18,'FontName','Arial')
box off
ylim([0 0.45])
saveas(gcf,[saveDir filesep 'Fig4B_MeanMag_vs_Time_DividingCell.png'])
saveas(gcf,[saveDir filesep 'Fig4B_MeanMag_vs_Time_DividingCell.svg'])

figure(2)
set(gcf,'Position',[400 500 560*1.5 420])
plot(t,meanmag(2,4:75),'Color',0.6*[1 1 1],'LineWidth',0.5)
hold on
scatter(t,meanmag(2,4:75),50,meanmag(2,4:75),'filled')
hold off
colormap(colorcet('L8'))
clim([0.15 0.4]) % This is different than the quiver graphs! The mean will be lower because there are always regions that are not moving.
ylabel(['Mean Magnitude (' char(181) 'm/min)'])
xlabel('Time (min)')
set(gca,'FontSize',18,'FontName','Arial')
box off
ylim([0 0.45])
saveas(gcf,[saveDir filesep 'Fig4E_MeanMag_vs_Time_MigratingCell.png'])
saveas(gcf,[saveDir filesep 'Fig4E_MeanMag_vs_Time_MigratingCell.svg'])


%% Phi Distribution - Dividing Cell

saturation = 0.03;

phiC = repmat(phiCenters(:),[1, length(t)]);
phiAlpha = distPhiW(:,4:75,1);
maxVal = max(phiAlpha(:));
phiAlpha = (phiAlpha-min(phiAlpha(:)))/(max(phiAlpha(:))-min(phiAlpha(:)));
phiAlpha = phiAlpha*maxVal/saturation; % Saturate the colors a bit for clarity
phiAlpha(phiAlpha>1) = 1;

figure(3)
set(gcf,'Position',[400 500 560*1.5 420])
h = imagesc(t,phiBins*180/pi,phiC);
colormap(colorcet('D2'))
ylabel('\phi (degrees)')
xlabel('Time (min)')
set(gca,'FontSize',18,'FontName','Arial')
set(gca,'YTick',-90:45:90)
set(gca,'YDir','normal')
set(h,'alphadata',phiAlpha)
set(gca,'Color','k')
set(gcf,'Color','w')
box off
set(gcf, 'InvertHardCopy', 'off');

saveas(gcf,[saveDir filesep 'Fig4C_WeightedPhiDist_vs_Time_DividingCell.png'])
saveas(gcf,[saveDir filesep 'Fig4C_WeightedPhiDist_vs_Time_DividingCell.svg'])

%% Phi Distribution - Migrating Cell

phiC = repmat(phiCenters(:),[1, length(t)]);
phiAlpha = distPhiW(:,4:75,2);
maxVal = max(phiAlpha(:));
phiAlpha = (phiAlpha-min(phiAlpha(:)))/(max(phiAlpha(:))-min(phiAlpha(:)));
phiAlpha = phiAlpha*maxVal/saturation; % Saturate the colors a bit for clarity
phiAlpha(phiAlpha>1) = 1;

figure(4)
set(gcf,'Position',[400 500 560*1.5 420])
h = imagesc(t,phiBins*180/pi,phiC);
colormap(colorcet('D2'))
ylabel('\phi (degrees)')
xlabel('Time (min)')
set(gca,'FontSize',18,'FontName','Arial')
set(gca,'YTick',-90:45:90)
set(gca,'YDir','normal')
set(h,'alphadata',phiAlpha)
set(gca,'Color','k')
set(gcf,'Color','w')
box off
set(gcf, 'InvertHardCopy', 'off');

saveas(gcf,[saveDir filesep 'Fig4F_WeightedPhiDist_vs_Time_MigratingCell.png'])
saveas(gcf,[saveDir filesep 'Fig4F_WeightedPhiDist_vs_Time_MigratingCell.svg'])

%% Phi Distribution - Color Bar

cAlpha = linspace(0,1,length(t));
cAlpha = repmat(cAlpha,[size(phiC,1) 1]);
cRange = linspace(0,saturation,length(t)); % This is based on the saturation chosen above

figure(5)
set(gcf,'Position',[1200 500 560/2.5 420])
h = imagesc(cRange,phiBins*180/pi,phiC);
set(h,'alphadata',cAlpha)
colormap(colorcet('D2'))
set(gca,'YDir','normal')
set(gca,'FontSize',18,'FontName','Arial')
set(gca,'YTick',-90:45:90)
set(gca,'XTick',0:.02:max(cRange))
set(gca,'Color','k')
set(gcf,'Color','w')
ylabel('\phi (degrees)')
% xlabel({'Weighted';'Probability'})
xlabel('Probability')
set(gcf, 'InvertHardCopy', 'off');

saveas(gcf,[saveDir filesep 'Fig4C_WeightedPhiDist_ColorScaleBar.png'])
saveas(gcf,[saveDir filesep 'Fig4C_WeightedPhiDist_ColorScaleBar.svg'])

