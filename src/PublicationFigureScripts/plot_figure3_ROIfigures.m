clc
clear
close all

disp(datetime('now'))

% Paths
imDir = 'X:\Force Project\PublicationData\AiryScan_twoChannels\';
saveDir = [imDir '\PublicationFigures\Figure3\'];
movieDir = [imDir '\PublicationFigures\MovieS4\'];
ch0Name = 'C2-04302024_U2OS_Cherry-Myo_Lifeact-GFP_30s_3slices_2-Airyscan Processing-13'; % myosin
ch1Name = 'C1-04302024_U2OS_Cherry-Myo_Lifeact-GFP_30s_3slices_2-Airyscan Processing-13'; % actin
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

%%% ROI: x1 x2 y1 y2
rois = [50 200 300 450 ; ...% protrusion at front
    375 525 25 175]; % protrusion at top
Nrois = size(rois,1);

if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
if ~exist(movieDir,'dir')
    mkdir(movieDir)
end

cmapROI = [0, 0, 0; ... black
    230, 159, 0; ... orange
    86, 180, 233; ...sky blue
    0, 158, 155; ...bluish green
    240, 228, 66; ...yellow
    0, 114, 178; ... blue
    213, 94, 0;... vermillion
    204, 121, 167]/255; %reddish purple


%% Get Image Sizes
relFile = [ch0Name '_rel_t' num2str(4,'%04u') '.tiff'];
meta = imfinfo([ch0Dir filesep relFile]);
Nx = meta(1).Width;
Nx = Nx-2*crop;
Ny = meta(2).Height;
Ny = Ny-2*crop;
Nz = length(meta);

ch0_theta = NaN*ones(max(tRange),Nrois+1);
ch1_theta = ch0_theta;
ch0_theta_weight = ch0_theta;
ch1_theta_weight = ch1_theta;

%% Loop through frames
for exampleFrame = tRange

    disp([char(datetime('now')) ' - Processing frame ' num2str(exampleFrame) '...'])

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

    % %%% Check ROIs
    % figure(1)
    % imshow(sum(relMask,3)>0)
    % hold on
    % for kk = 1:Nrois
    %     plot([rois(kk,1) rois(kk,1) rois(kk,2) rois(kk,2) rois(kk,1) rois(kk,1)],...
    %         [rois(kk,3) rois(kk,4) rois(kk,4) rois(kk,3) rois(kk,3) rois(kk,4)],...
    %         'LineWidth',4)
    % end
    % hold off
    % % error('stop')


    %% Load in velocities
    
    vx = TIFFvolume([ch0Dir filesep ch0Name '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx = vx(crop+1:end-crop,crop+1:end-crop,:);
    vx = vx.*relMask./relMask*xyscale/tscale;
    
    vy = TIFFvolume([ch0Dir filesep ch0Name '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy = vy(crop+1:end-crop,crop+1:end-crop,:);
    vy = vy.*relMask./relMask*xyscale/tscale;

    % vz = TIFFvolume([ch0Dir filesep ch0Name '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    % vz = vz(crop+1:end-crop,crop+1:end-crop,:);
    % vz = vz.*relMask./relMask*zscale/tscale;

    vx1 = TIFFvolume([ch1Dir filesep ch1Name '_vx_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vx1 = vx1(crop+1:end-crop,crop+1:end-crop,:);
    vx1 = vx1.*relMask./relMask*xyscale/tscale;
    
    vy1 = TIFFvolume([ch1Dir filesep ch1Name '_vy_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    vy1 = vy1(crop+1:end-crop,crop+1:end-crop,:);
    vy1 = vy1.*relMask./relMask*xyscale/tscale;
    % 
    % vz1 = TIFFvolume([ch1Dir filesep ch1Name '_vz_t' num2str(exampleFrame,'%04u') '.tiff'],Nz);
    % vz1 = vz1(crop+1:end-crop,crop+1:end-crop,:);
    % vz1 = vz1.*relMask./relMask*zscale/tscale;

    theta = atan2(vy,vx);   
    theta1 = atan2(vy1,vx1);

    % mag = sqrt(vx.^2 + vy.^2 + vz.^2);
    % mag1 = sqrt(vx1.^2 + vy1.^2 + vz1.^2);
    mag = sqrt(vx.^2 + vy.^2);
    mag1 = sqrt(vx1.^2 + vy1.^2);

    %% Means
    badRows = isnan(theta) | isnan(theta1);

    term1 = mean(sin(theta(~badRows)));
    term2 = mean(cos(theta(~badRows)));
    ch0_theta(exampleFrame,1) = atan2(term1,term2);

    term1 = mean(sin(theta1(~badRows)));
    term2 = mean(cos(theta1(~badRows)));
    ch1_theta(exampleFrame,1) = atan2(term1,term2);

    term1 = mean(mag(~badRows).*sin(theta(~badRows)))/sum(mag(~badRows));
    term2 = mean(mag(~badRows).*cos(theta(~badRows)))/sum(mag(~badRows));
    ch0_theta_weight(exampleFrame,1) = atan2(term1,term2);

    term1 = mean(mag1(~badRows).*sin(theta1(~badRows)))/sum(mag1(~badRows));
    term2 = mean(mag1(~badRows).*cos(theta1(~badRows)))/sum(mag1(~badRows));
    ch1_theta_weight(exampleFrame,1) = atan2(term1,term2);

    for jj = 1:Nrois

        thetaSlice = theta(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        thetaSlice1 = theta1(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        badSlice = isnan(thetaSlice) | isnan(thetaSlice1);

        slice=theta(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        sliceMag = mag(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        term1 = mean(sin(slice(~badSlice)));
        term2 = mean(cos(slice(~badSlice)));
        ch0_theta(exampleFrame,jj+1) = atan2(term1,term2);

        term1 = mean(sliceMag(~badSlice).*sin(slice(~badSlice)))/sum(sliceMag(~badSlice));
        term2 = mean(sliceMag(~badSlice).*cos(slice(~badSlice)))/sum(sliceMag(~badSlice));
        ch0_theta_weight(exampleFrame,jj+1) = atan2(term1,term2);

        slice=theta1(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        sliceMag = mag1(rois(jj,3):rois(jj,4),rois(jj,1):rois(jj,2),:);
        term1 = mean(sin(slice(~badSlice)));
        term2 = mean(cos(slice(~badSlice)));
        ch1_theta(exampleFrame,jj+1) = atan2(term1,term2);

        term1 = mean(sliceMag(~badSlice).*sin(slice(~badSlice)))/sum(sliceMag(~badSlice));
        term2 = mean(sliceMag(~badSlice).*cos(slice(~badSlice)))/sum(sliceMag(~badSlice));
        ch1_theta_weight(exampleFrame,jj+1) = atan2(term1,term2);

    end  

    %% Clean up
    % close all
    disp([char(datetime('now')) ' - Frame ' num2str(exampleFrame) ' complete'])

end

%% Overall Mean Theta

tVec = (tRange-1)*tscale;
tVec2 = (0:(max(tRange)+2))*tscale;

% Using negative theta for the convention that up = +90 and down = -90
figure(2)
set(gcf,'Position',[550 500 500 500])
plot(tVec,-ch0_theta_weight(tRange,1)*180/pi,'-o','Color',0.3*[1 1 1],'LineWidth',2,'MarkerFaceColor',0.3*[1 1 1],'MarkerSize',10)
hold on
plot(tVec,-ch1_theta_weight(tRange,1)*180/pi,'-s','Color',0.7*[1 1 1],'LineWidth',2,'MarkerFaceColor',0.7*[1 1 1],'MarkerSize',10)
plot(tVec2,0*tVec2,'--k')
hold off
xlabel('Time (min)')
ylabel('Mean \theta')
set(gca,'FontSize',20)
box off
legend('Myosin','Actin','Location','NorthWest','Orientation','Horizontal')
legend boxoff
ylim([-180 180])
set(gca,'YTick',-180:90:180)

saveas(gcf,[saveDir filesep 'Fig3E_MeanTheta.png'])
saveas(gcf,[saveDir filesep 'Fig3E_MeanTheta.svg'])


%% Theta in ROIs

for kk = 1:Nrois

    figure(kk+2)
    set(gcf,'Position',[550 500 500 400])

    plot(tVec,-180/pi*ch0_theta_weight(tRange,kk+1),'-o','Color',0.6*cmapROI(kk+1,:),'MarkerFaceColor',0.6*cmapROI(kk+1,:),'MarkerSize',10,'LineWidth',2);
    hold on
    plot(tVec,-180/pi*ch1_theta_weight(tRange,kk+1),'-s','Color',cmapROI(kk+1,:),'MarkerFaceColor',cmapROI(kk+1,:),'MarkerSize',10,'LineWidth',2);
    plot(tVec2,0*tVec2,'--k')
    hold off
    xlabel('Time (min)')
    ylabel('Mean \theta')
    set(gca,'FontSize',20)
    box off
    legend('Myosin','Actin','Location','NorthWest','Orientation','Horizontal')
    legend boxoff
    ylim([-180 180])
    set(gca,'YTick',-180:90:180)
    
    saveas(gcf,[saveDir filesep 'Fig3F-G_MeanTheta_ROI' num2str(kk) '.png'])
    saveas(gcf,[saveDir filesep 'Fig3F-G_MeanTheta_ROI' num2str(kk) '.svg'])

end


%% Image Panels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at the photobleaching

myo = TIFFvolume([imDir filesep ch0Name '.tif'],31*3);
myo = reshape(myo,[size(myo,1), size(myo,2), 3, 31]);
myo = myo(crop+1:end-crop,crop+1:end-crop,:,:);

act = TIFFvolume([imDir filesep ch1Name '.tif'],31*3);
act = reshape(act,[size(act,1), size(act,2), 3, 31]);
act = act(crop+1:end-crop,crop+1:end-crop,:,:);


myoM = NaN*ones(size(myo,4),1);
actM = myoM;

for kk = 1:size(myo,4)
    slice = myo(:,:,:,kk);
    myoM(kk) = mean(slice(~isnan(slice)));

    slice = act(:,:,:,kk);
    actM(kk) = mean(slice(~isnan(slice)));
end

myoDF = (myoM-myoM(1))/myoM(1);
actDF = (actM-actM(1))/actM(1);

figure(5)
plot((0:length(myoM)-1)*tscale,myoDF*100,'Color',[0.1 0.9 0.1],'LineWidth',3)
hold on
plot((0:length(myoM)-1)*tscale,actDF*100,'Color',[0.9 0.1 0.9],'LineWidth',3)
plot((0:length(myoM)-1)*tscale,0*(1:length(myoM))*tscale,'--k')
hold off
set(gca,'FontSize',20)
box off
xlabel('Time (min)')
ylabel('\DeltaF/F (%)')
legend('Myosin','Actin','Location','SouthWest')
legend boxoff

saveas(gcf,[saveDir filesep 'Fig3Supp_dfF_vs_Time.png'])
saveas(gcf,[saveDir filesep 'Fig3Supp_dfF_vs_Time.svg'])

%% ROI Images for Figure

A = squeeze(max(myo,[],3));
imNow = A(:,:,1);
imNow = imadjust(imNow,[25 4000]/2^16,[],0.4);

for kk = 1:Nrois
    imNow = insertShape(imNow,'Rectangle',[rois(kk,1) rois(kk,3) rois(kk,2)-rois(kk,1) rois(kk,4)-rois(kk,3)],'LineWidth',6,'Color',cmapROI(kk+1,:)*2^16-1);
end
a = size(imNow,1)-30;
imNow = insertShape(imNow,'Line',[20 a 20+5/xyscale a],'Color',2^16*[1 1 1],'LineWidth',20);

figure(6)
imshow(imNow)

imwrite(imNow,[saveDir filesep 'Fig3D_myosin_with_ROIs_scaleBar5um.png'])

B = squeeze(max(act,[],3));
imNow = B(:,:,1);
imNow = imadjust(imNow,[25 4000]/2^16);

for kk = 1:Nrois
    imNow = insertShape(imNow,'Rectangle',[rois(kk,1) rois(kk,3) rois(kk,2)-rois(kk,1) rois(kk,4)-rois(kk,3)],'LineWidth',6,'Color',cmapROI(kk+1,:)*2^16-1);
end
imNow = insertShape(imNow,'Line',[20 a 20+5/xyscale a],'Color',2^16*[1 1 1],'LineWidth',20);

figure(7)
imshow(imNow)

imwrite(imNow,[saveDir filesep 'Fig3D_actin_with_ROIs_scaleBar5um.png'])

%% ROI Snapshots
% Use two frames to illustrate the ROIs

t1 = 16;
t2 = 17;

disp(['t1 = ' num2str(t1) ' = ' num2str(t1*tscale) ' min'])
disp(['t2 = ' num2str(t2) ' = ' num2str(t2*tscale) ' min'])

for kk = 1 % Only look at the ROI with opposing motion

    %%%% Myosin

    roi1 = A(rois(kk,3):rois(kk,4),rois(kk,1):rois(kk,2),t1);
    roi1 = imadjust(roi1,[25 4000]/2^16,[],0.4);

    roi2 = A(rois(kk,3):rois(kk,4),rois(kk,1):rois(kk,2),t2);
    roi2 = imadjust(roi2,[25 4000]/2^16,[],0.4);

    clear roi3
    roi3(:,:,1) = roi2;
    roi3(:,:,2) = roi1;
    roi3(:,:,3) = roi2;

    imwrite(roi3,[saveDir filesep 'Fig3H_ROI' num2str(kk) '_myosin_t1green_t2magneta.png'])

    %%%% Actin

    roi1 = B(rois(kk,3):rois(kk,4),rois(kk,1):rois(kk,2),t1);
    roi1 = imadjust(roi1,[25 4000]/2^16);

    roi2 = B(rois(kk,3):rois(kk,4),rois(kk,1):rois(kk,2),t2);
    roi2 = imadjust(roi2,[25 4000]/2^16);

    clear roi3
    roi3(:,:,1) = roi2;
    roi3(:,:,2) = roi1;
    roi3(:,:,3) = roi2;

    imwrite(roi3,[saveDir filesep 'Fig3H_ROI' num2str(kk) '_actin_t1green_t2magneta.png'])

end


%% Clean up
close all

%% Corresponding Movie S4

for exampleFrame = 1:size(A,3)

    imNowA = A(:,:,exampleFrame);
    imNowA = imadjust(imNowA,[25 4000]/2^16,[],0.4);
    imNowA = uint8(double(imNowA)/(2^16-1)*255);
    
    imNowB = B(:,:,exampleFrame);
    imNowB = imadjust(imNowB,[25 4000]/2^16);
    imNowB = uint8(double(imNowB)/(2^16-1)*255);
    
    imNowC(:,:,1) = imNowB;
    imNowC(:,:,2) = imNowA;
    imNowC(:,:,3) = imNowB;
    
    for kk = 1:Nrois
        imNowA = insertShape(imNowA,'Rectangle',[rois(kk,1) rois(kk,3) rois(kk,2)-rois(kk,1) rois(kk,4)-rois(kk,3)],'LineWidth',6,'Color',cmapROI(kk+1,:)*255);
    end

    % Add annotations
    a = size(imNowA,1)-30;
    imNowA = insertShape(imNowA,'Line',[20 a 20+5/xyscale a],'Color',2^16*[1 1 1],'LineWidth',20); % scalebar
    tText = (exampleFrame-1)*tscale;
    mm = floor(tText);
    ss = tText-mm;
    ss = round(ss*60); % seconds
    tText = [num2str(mm,'%02u') ':' num2str(ss,'%02u')];
    imNowA = insertText(imNowA,[20 15],tText,BoxOpacity=0,TextColor='white',FontSize=50);
    
    for kk = 1:Nrois
        imNowB = insertShape(imNowB,'Rectangle',[rois(kk,1) rois(kk,3) rois(kk,2)-rois(kk,1) rois(kk,4)-rois(kk,3)],'LineWidth',6,'Color',cmapROI(kk+1,:)*255);
    end
    
    for kk = 1:Nrois
        imNowC = insertShape(imNowC,'Rectangle',[rois(kk,1) rois(kk,3) rois(kk,2)-rois(kk,1) rois(kk,4)-rois(kk,3)],'LineWidth',6,'Color',cmapROI(kk+1,:)*255);
    end
    
    figure(1)
    imshow([imNowA imNowB imNowC])
    
    imwrite([imNowA imNowB imNowC],[movieDir filesep 'threePanels_frame' num2str(exampleFrame,'%03u') '_scaleBar5um.png'])

end

%% Clean up
close all
disp('Movie Complete')