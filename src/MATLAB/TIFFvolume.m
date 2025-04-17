function vol = TIFFvolume(filename,Nx,Ny,frames)
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
%
% OUTPUTS:
% vol       = Matrix of size [Ny,Nx,frames] and type double
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vol = ones(Ny,Nx,frames);
    intiff = Tiff(filename,"r");
    for ii = 1:frames-1
        vol(:,:,ii) = read(intiff);
        nextDirectory(intiff);
    end
    vol(:,:,frames) = read(intiff);

end