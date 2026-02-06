% Start with a clean workspace
clc, clear, close all
profile -memory on

%%% Set up optical flow parameters
xyzSig = 3; % Spatial smoothing (voxels)
tSig = 1; % Temporal smoothing (frames)
wSig = 4; % Lucas-Kanade neighborhood size (voxels)

% Paths
mainDir = 'X:\Force Project\PublicationData\PerformanceTests';
saveName = [mainDir filesep 'PerformanceData_MATLAB_Laptop.csv'];
program = {'MATLAB-laptop'};


%% 2D OneTif

% Set up
imDir = [mainDir filesep '2D_OneTif'];
fileType = 'OneTif';
spatialDimensions = 2;

% Get the list of tifs in this folder
tifList = dir([imDir filesep '*.tif']);

profile clear

for jj = 1:length(tifList) % Loop through all tifs

    imName = tifList(jj).name(1:end-4);
    disp(['FILE: ' imName])

    meta = imfinfo([imDir filesep imName '.tif'],'tif');
    Nx = meta(1).Width;
    Ny = meta(1).Height;
    Nvoxels = Nx*Ny;
    clear Nx Ny meta

    for ii = 1:10 % Iterations to get statistics on the run time and memory

        disp(['Iteration ' num2str(ii)])
    
        % 
        profile on
        
        % Run the optical flow
        process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)
        
        p = profile('info');
    
        % Find the upper-level function (all we care about for this particular anlaysis)
        for mm = 1:length(p.FunctionTable)
            if strcmp('process_flow',p.FunctionTable(mm).FunctionName)
                break
            end
        end
        
        % Save the results
        TotalTime = p.FunctionTable(mm).TotalTime;
        PeakMem = p.FunctionTable(mm).PeakMem;
        imageDirectory = {imDir};
        imageName = {imName};
        processingType = {fileType};
        param = table(program,imageDirectory,imageName,processingType,spatialDimensions,xyzSig,tSig,wSig,TotalTime,PeakMem,Nvoxels);
        writetable(param,saveName,'WriteMode','Append')
        
        % Clear before next round
        profile clear
    
    end

end

clear tifList

%% 2D Sequence

% Set up
imDir = [mainDir filesep '2D_Sequence'];
imName = '20241107_U2OS_SGRLC_100Xoil_15mintimelapse_01_25plaser-slice7_t*';
fileType = 'SequenceT';
spatialDimensions = 2;

disp(['FILE: ' imName])

tifList = dir([imDir filesep imName '.tif']);
meta = imfinfo([imDir filesep tifList(1).name],'tif');
Nx = meta(1).Width;
Ny = meta(1).Height;
Nvoxels = Nx*Ny;
clear Nx Ny meta

for ii = 1:10 % Iterations to get statistics on the run time and memory

    disp(['Iteration ' num2str(ii)])

    % profile -memory on
    profile on
    
    % Run the optical flow
    process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)
    
    p = profile('info');

    % Find the upper-level function (all we care about for this particular anlaysis)
    for mm = 1:length(p.FunctionTable)
        if strcmp('process_flow',p.FunctionTable(mm).FunctionName)
            break
        end
    end
    
    % Save the results
    TotalTime = p.FunctionTable(mm).TotalTime;
    PeakMem = p.FunctionTable(mm).PeakMem;
    imageDirectory = {imDir};
    imageName = {imName};
    processingType = {fileType};
    param = table(program,imageDirectory,imageName,processingType,spatialDimensions,xyzSig,tSig,wSig,TotalTime,PeakMem,Nvoxels);
    writetable(param,saveName,'WriteMode','Append')
    
    % Clear before next round
    profile clear

end


%% 3D OneTif

% Set up
imDir = [mainDir filesep '3D_OneTif'];
fileType = 'OneTif';
spatialDimensions = 3;

% Get the list of tifs in this folder
tifList = dir([imDir filesep '*.tif']);

for jj = 2:length(tifList) % Skip the larger z-stack for the laptop tests

    imName = tifList(jj).name(1:end-4);
    disp(['FILE: ' imName])

    meta = imfinfo([imDir filesep imName '.tif'],'tif');
    Nx = meta(1).Width;
    Ny = meta(1).Height;
    Nz = meta(1).ImageDescription;
    Nz = strsplit(Nz,'\n');
    Nz = strsplit(Nz{3},'=');
    Nz = str2double(Nz{end});
    Nvoxels = Nx*Ny*Nz;
    clear Nx Ny Nz meta

    for ii = 1:10 % Iterations to get statistics on the run time and memory

        disp(['Iteration ' num2str(ii)])
    
        % profile -memory on
        profile on
        
        % Run the optical flow
        process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)
        
        p = profile('info');
    
        % Find the upper-level function (all we care about for this particular anlaysis)
        for mm = 1:length(p.FunctionTable)
            if strcmp('process_flow',p.FunctionTable(mm).FunctionName)
                break
            end
        end
        
        % Save the results
        TotalTime = p.FunctionTable(mm).TotalTime;
        PeakMem = p.FunctionTable(mm).PeakMem;
        imageDirectory = {imDir};
        imageName = {imName};
        processingType = {fileType};
        param = table(program,imageDirectory,imageName,processingType,spatialDimensions,xyzSig,tSig,wSig,TotalTime,PeakMem,Nvoxels);
        writetable(param,saveName,'WriteMode','Append')
        
        % Clear before next round
        profile clear
    
    end

end

%% 3D Sequence
% 
% Set up
% imDir = [mainDir filesep '3D_Sequence'];
% imName = 'scan_Cam1_ch0_tile0_t*_deskew_after_decon';
% fileType = 'SequenceT';
% spatialDimensions = 3;
% 
% disp(['FILE: ' imName])
% 
% tifList = dir([imDir filesep imName '.tif']);
% meta = imfinfo([imDir filesep tifList(1).name],'tif');
% Nx = meta(1).Width;
% Ny = meta(1).Height;
% Nz = length(meta);
% Nvoxels = Nx*Ny*Nz;
% clear Nx Ny Nz meta
% 
% for ii = 1:10 % Iterations to get statistics on the run time and memory
% 
%     disp(['Iteration ' num2str(ii)])
% 
%     profile -memory on
%     profile on
% 
%     Run the optical flow
%     process_flow(imDir,imName,fileType,spatialDimensions,xyzSig,tSig,wSig)
% 
%     p = profile('info');
% 
%     Find the upper-level function (all we care about for this particular anlaysis)
%     for mm = 1:length(p.FunctionTable)
%         if strcmp('process_flow',p.FunctionTable(mm).FunctionName)
%             break
%         end
%     end
% 
%     Save the results
%     TotalTime = p.FunctionTable(mm).TotalTime;
%     PeakMem = p.FunctionTable(mm).PeakMem;
%     imageDirectory = {imDir};
%     imageName = {imName};
%     processingType = {fileType};
%     param = table(program,imageDirectory,imageName,processingType,spatialDimensions,xyzSig,tSig,wSig,TotalTime,PeakMem,Nvoxels);
%     writetable(param,saveName,'WriteMode','Append')
% 
%     Clear before next round
%     profile clear
% 
% end

%%
clc
disp('Processing Complete')