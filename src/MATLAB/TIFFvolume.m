% Use TIFF to load in 3D volumes
function vol = TIFFvolume(filename,Nx,Ny,frames)

    vol = ones(Ny,Nx,frames);
    intiff = Tiff(filename,"r");
    for ii = 1:frames-1
        vol(:,:,ii) = read(intiff);
        nextDirectory(intiff);
    end
    vol(:,:,frames) = read(intiff);

end