% This script runs the optical flow over the data presented in the
% publication. It does not create visualizations or run any analyses.

% Start with a clean workspace
clc, clear, close all

%%% Set up optical flow parameters
xyzSig = 3; % Spatial smoothing (voxels)
tSig = 1; % Temporal smoothing (frames)
wSig = 4; % Lucas-Kanade neighborhood size (voxels)

% Paths
imDir = 'X:\Force Project\PublicationData\MOSAIC_Actin\deskew_after_decon'; 
imName = 'scan_Cam1_ch0_tile0_t*_deskew_after_decon';
savedir = 'X:\Force Project\PublicationData\MOSAIC_Actin\PublicationFigures\FigureS2';

% Metadata
xyscale = 0.108; % um/pixel
zscale = 0.5*sin(32.45*pi/180); % um/pixel
tscale = 60075/1000/60; % minutes/frame
Nt = 78;

% For Processing
dt = 1:10; % Range of dt to consider
% Crop to a region around the dividing cell to avoid unncessary processing
% time for this illustration
y1 = 500;
y2 = 1100;
x1 = 900;
x2 = 1500;
z1 = 120;
z2 = 300;
% Size after cropping
Nx = 601;
Nz = 181;
% Frame to use
Ntslice = floor(Nt/2)-1;

% Visualization
zSlice = 43;
gap = 6;
qscale = 4;

%% Run the flow, looping through dt

list = dir([imDir filesep imName '.tif']);

for kk = 1:length(dt)

    loopStart = datetime('now');
    disp([char(datetime('now')) ' - Processing dt ' num2str(dt(kk)) '...'])

    dtNow = dt(kk);
    tNow = Ntslice+dtNow*(-3:3); % This keeps the center frame always the same, but changes dt

    savenow = [savedir filesep num2str(dtNow)];
    if ~exist(savenow,'dir')
        mkdir(savenow)
    end

    images = NaN*ones(Nx,Nx,Nz,7);
    for jj = 1:length(tNow)

        tmp = TIFFvolume([imDir filesep list(tNow(jj)).name],z2);
        images(:,:,:,jj) = tmp(y1:y2,x1:x2,z1:end);
        clear tmp

    end

    images = double(images); % For calculations

    %  Run the optical flow
    [vx,vy,vz,rel] = calc_flow3D(images ,xyzSig, tSig, wSig);
    clear images

    % Save this frame
    TIFFwrite([savenow filesep 'dt' num2str(dtNow) '_vx_t' num2str(Ntslice,'%04u') '.tiff'],vx)
    TIFFwrite([savenow filesep 'dt' num2str(dtNow) '_vy_t' num2str(Ntslice,'%04u') '.tiff'],vy)
    TIFFwrite([savenow filesep 'dt' num2str(dtNow) '_vz_t' num2str(Ntslice,'%04u') '.tiff'],vz)
    TIFFwrite([savenow filesep 'dt' num2str(dtNow) '_rel_t' num2str(Ntslice,'%04u') '.tiff'],rel)
    clear rel vx vy vz

    framestime = datetime('now');
    disp([char(datetime('now')) ' - dt ' num2str(dt(kk)) ' saved.  Duration: ' char(framestime-loopStart)])

end

%% Make quiver figures

x = (0:Nx-1)*xyscale;
z = (0:Nz-1)*zscale;
[X,Y,Z] = meshgrid(x,x,z);

sMax = 0.75;
bins = linspace(0,sMax,255);
cmap = colorcet('L8');
cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

xNow = X(1:gap:end,1:gap:end,zSlice);
yNow = Y(1:gap:end,1:gap:end,zSlice);

im1 = TIFFvolume([imDir filesep list(Ntslice).name],z2);
im1 = im1(y1:y2,x1:x2,z1:end);

for kk = 1:length(dt)

    tic

    dtNow = dt(kk);
    disp(dtNow)

    vx = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vx_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    vy = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vy_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    vz = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vz_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    rel = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_rel_t' num2str(Ntslice,'%04u') '.tiff'],Nz);

    relThresh = prctile(double(rel(:)),92);
    relMask = rel > relThresh;

    im2 = TIFFvolume([imDir filesep list(Ntslice+dtNow).name],z2);
    im2 = im2(y1:y2,x1:x2,z1:end);

    vx = vx.*relMask*xyscale/tscale;
    vy = vy.*relMask*xyscale/tscale;
    vz = vz.*relMask*xyscale/tscale;

    speed = sqrt(vx.^2 + vy.^2 + vz.^2);  

    figure(1)
    C = imfuse(imadjust(im1(:,:,zSlice),[0 0.006]),imadjust(im2(:,:,zSlice),[0 0.006]));
    C = insertShape(C,'Line',[20, 570, 20+10/xyscale 570],'Color','w','LineWidth',20);
    imshow(C)
    imwrite(C,[savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_ColorMerge.png'])    
    
    vxNow = vx(1:gap:end,1:gap:end,zSlice);
    vyNow = vy(1:gap:end,1:gap:end,zSlice);
    speedNow = speed(1:gap:end,1:gap:end,zSlice);
    
    figure(2)
    set(gcf,'Position',[500 250 1100 900])
    for jj = 1:length(bins)-1     
        slice = (speedNow>=bins(jj)) & (speedNow<bins(jj+1));
        j = quiver(xNow(slice),yNow(slice),vxNow(slice),vyNow(slice),0,...
            'Color',cmap(jj,:),'LineWidth',1);
        hold on
        if qscale ~= 0
            hU = get(j,'UData') ;
            hV = get(j,'VData') ;
            set(j,'UData',qscale*hU,'VData',qscale*hV)
        end
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',20)
    set(gca,'XTickLabel',{},'YTickLabel',{})
    set(gca,'ydir','reverse')
    set(gca,'Color','k')
    xlim([min(xNow(:)) max(xNow(:))])
    ylim([min(yNow(:)) max(yNow(:))])
    colormap(cmap)
    h = colorbar('Location','SouthOutside');
    set(h,'Ticks',0:.25:sMax)
    set(get(h,'Label'),'String',['Magnitude (' char(181) 'm/min)'])
    clim([0 sMax])
    set(gca,'Color','k')
    set(gcf,'Color','w')
    set(gcf, 'InvertHardCopy', 'off');
    
    f = getframe(gca);
    imwrite(f.cdata,[savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_Magnitude_3D_zSlice' num2str(zSlice) '.png'])
    saveas(gcf,[savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_Magnitude_3D_zSlice' num2str(zSlice) '_withColorbar.png'])

    toc

end

%% Calculate Distributions 

binsV = linspace(-3,3,255);
binsS = linspace(0,0.5,255);

Nvx = NaN*ones(length(dt),length(binsV)-1);
Nvy = Nvx;
Nvz = Nvx;

Ns = NaN*ones(length(dt),length(binsS)-1);
Nsdt = Ns;

tic
for kk = 1:length(dt)

    dtNow = dt(kk);
    disp(dtNow)

    vx = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vx_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    vy = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vy_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    vz = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_vz_t' num2str(Ntslice,'%04u') '.tiff'],Nz);
    rel = TIFFvolume([savedir filesep num2str(dtNow) filesep 'dt' num2str(dtNow) '_rel_t' num2str(Ntslice,'%04u') '.tiff'],Nz);

    relThresh = prctile(double(rel(:)),92);
    relMask = rel > relThresh;
    relMask = relMask./relMask; % Introduce NaN

    % Distributions in pixels
    vx = vx.*relMask;
    vy = vy.*relMask;
    vz = vz.*relMask;

    Nvx(kk,:) = histcounts(vx(:),binsV,'Normalization','Probability');
    Nvy(kk,:) = histcounts(vy(:),binsV,'Normalization','Probability');
    Nvz(kk,:) = histcounts(vz(:),binsV,'Normalization','Probability');

    % Speed in real units
    vx = vx*xyscale/tscale;
    vy = vy*xyscale/tscale;
    vz = vz*xyscale/tscale;

    speed = sqrt(vx.^2 + vy.^2 + vz.^2);
    % theta = atan2(-vy,vx);
    % phi = atan(vz./sqrt(vx.^2+vy.^2));

    Ns(kk,:) = histcounts(speed(:),binsS,'Normalization','Probability');
    Nsdt(kk,:) = histcounts(speed(:)/dtNow,binsS,'Normalization','Probability');

end
toc

%% Figure S2 E - Magnitude Distributions

figure(3);
set(gcf,'Position',[600 500 560 525])

t = tiledlayout('vertical');
cmap = colorcet('L19');
cmap = cmap(round(linspace(1,length(cmap),length(dt)+6)),:);
cmap = cmap(6:end,:);

nexttile
for kk = 1:length(dt)
    stairs(binsS,[Ns(kk,:), Ns(kk,end)],'LineWidth',2,'Color',cmap(kk,:))
    hold on
end
hold off
set(gca,'FontSize',20)
box off
xlim([min(binsS),max(binsS)])
set(gca,'XTick',0:.1:.5)
ylim([0 10^-3])
xlabel(['Magnitude (' char(181) 'm/min)'],'FontSize',20)

nexttile
for kk = 1:length(dt)
    stairs(binsS,[Nsdt(kk,:), Nsdt(kk,end)],'LineWidth',2,'Color',cmap(kk,:))
    hold on
end
set(gca,'FontSize',20)
box off
xlim([min(binsS),max(binsS)])
set(gca,'XTick',0:.1:.5)
ylim([0 4*10^-3])
xlabel('Magnitude /\Deltat','FontSize',20)

ylabel(t, 'Probability Density','FontSize',20)

saveas(gcf,[savedir filesep 'FigureS2E_MagnitudeDistributions.png'])
saveas(gcf,[savedir filesep 'FigureS2E_MagnitudeDistributions.svg'])

%% Figure S2 D - Velocity Distributions

figure(6);
set(gcf,'Position',[600 500 560 750])

t = tiledlayout('vertical');
cmap = colorcet('L19');
cmap = cmap(round(linspace(1,length(cmap),length(dt)+6)),:);
cmap = cmap(6:end,:);

nexttile
for kk = 1:length(dt)
    stairs(binsV,[Nvx(kk,:), Nvx(kk,end)],'LineWidth',2,'Color',cmap(kk,:))
    hold on
end
hold off
set(gca,'FontSize',20)
box off
xlim([min(binsV),max(binsV)])
title('x')
ylim([0, 1.4*10^-3])
set(gca,'XTick',-3:3)


nexttile
for kk = 1:length(dt)
    stairs(binsV,[Nvy(kk,:), Nvy(kk,end)],'LineWidth',2,'Color',cmap(kk,:))
    hold on
end
hold off
set(gca,'FontSize',20)
box off
xlim([min(binsV),max(binsV)])
title('y')
ylim([0, 1.4*10^-3])
set(gca,'XTick',-3:3)

nexttile
for kk = 1:length(dt)
    stairs(binsV,[Nvz(kk,:), Nvz(kk,end)],'LineWidth',2,'Color',cmap(kk,:))
    hold on
end
hold off
set(gca,'FontSize',20)
box off
xlim([min(binsV),max(binsV)])
title('z')
ylim([0, 1.4*10^-3])

% title(t, 'Common title')
xlabel(t, 'Magnitude (pixels)','FontSize',20)
ylabel(t, 'Probability Density','FontSize',20)
set(gca,'XTick',-3:3)

saveas(gcf,[savedir filesep 'FigureS2D_VelocityDistributions.png'])
saveas(gcf,[savedir filesep 'FigureS2D_VelocityDistributions.svg'])

%% Vector Graphic Color Bar

sMax = 0.75;
bins = linspace(0,sMax,255);
cmap = colorcet('L8');
cmap = cmap(round(linspace(1,length(cmap),length(bins))),:);

figure(7)
set(gcf,'Position',[600 500 560 800])
imagesc(rand(100))

set(gca,'FontSize',20)
set(gca,'XTickLabel',{},'YTickLabel',{})
set(gca,'ydir','reverse')
set(gca,'Color','k')
colormap(cmap)
h = colorbar('Location','SouthOutside');
set(h,'Ticks',0:.25:sMax)
set(get(h,'Label'),'String',['Magnitude (' char(181) 'm/min)'])
clim([0 sMax])
set(gca,'Color','k')
set(gcf,'Color','w')
set(gcf, 'InvertHardCopy', 'off');

saveas(gcf,[savedir filesep 'FigureS2A-C_Colorbar.svg'])
