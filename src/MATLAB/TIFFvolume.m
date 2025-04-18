function vol = TIFFvolume(filename,Nx,Ny,frames,startFrame)
% Use TIFF to load in 3D volumes
%
% TIFF files are read using the TIFF library.
%
% vol = TIFFvolume(filename,Nx,Ny,frames)
%
% INPUTS:
% filename  = path for the saved file
% Nx        = Size of the image in X (width)
% Ny        = Size of the image in Y (height)
% frames    = Number of frames from the volume to load. This can be time
%               and/or z depending on the input file.
% startFrame = First frame to load from. Allows for only loading part of
%               the image. Frames startFrame:startFrame+frames are loaded.
%               Default value = 1 (the first frame).
%
% OUTPUTS:
% vol       = Matrix of size [Ny,Nx,frames] and type double
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use default inputs.
if ~exist('startFrame','var') || isempty(startFrame)
    startFrame = 1;
end
if startFrame < 1
    error('TIFFvolume requires the startFrame >=1 (this function uses 1-based indexing)')
end

%  Open the file for reading
intiff = Tiff(filename,"r");

% If not starting from 1, progress to the correct directory
if startFrame > 1
    for ii = 1:startFrame-1
        nextDirectory(intiff);
    end
end

% Load the number of frames requested
vol = ones(Ny,Nx,frames);
for ii = 1:frames-1
    vol(:,:,ii) = read(intiff);
    nextDirectory(intiff);
end
vol(:,:,frames) = read(intiff);

% Close the file
intiff.close();