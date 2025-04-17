% Use TIFF to save float tifs
function TIFFwrite(filename,A)

    % TIFF Tags
    tagstruct.ImageLength = size(A,1);
    tagstruct.ImageWidth = size(A,2);
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.BitsPerSample = 64;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Compression = Tiff.Compression.LZW;

    outtiff = Tiff(filename,'w8'); % Using bigtiff writing
    for ii = 1:size(A,3)-1
        outtiff.setTag(tagstruct);
        outtiff.write(A(:,:,ii));
        outtiff.writeDirectory();
    end
    outtiff.setTag(tagstruct);
    outtiff.write(A(:,:,size(A,3)));

    % Close the file
    outtiff.close();
    
end