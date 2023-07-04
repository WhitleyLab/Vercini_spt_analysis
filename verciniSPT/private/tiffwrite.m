function tiffwrite(filename,timg)
%function tiffwrite(filename,timg)

%decide the image type
imclass = class(timg);
if strcmp(imclass,'single')
    BitsPerSample = 32; 
    SampleFormat = Tiff.SampleFormat.IEEEFP; 
elseif strcmp(imclass,'double')
    timg = cast(timg,'single');
    BitsPerSample = 32; 
    SampleFormat = Tiff.SampleFormat.IEEEFP; 
elseif strcmp(imclass,'uint16')
    BitsPerSample = 16; 
    SampleFormat = Tiff.SampleFormat.UInt; 
elseif strcmp(imclass,'uint8')
    BitsPerSample = 8; 
    SampleFormat = Tiff.SampleFormat.UInt; 
else
    error('Unsupported image type');
end

t = Tiff(filename, 'w'); 
tagstruct.ImageLength = size(timg, 1); 
tagstruct.ImageWidth = size(timg, 2); 
nframe = size(timg, 3); 

tagstruct.Compression = Tiff.Compression.None; 
tagstruct.SampleFormat = SampleFormat;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
tagstruct.BitsPerSample = BitsPerSample;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
t.setTag(tagstruct); 
t.write(timg(:,:,1)); 

if nframe > 1
   for ii=2:nframe
        t.writeDirectory();
        t.setTag(tagstruct);
        t.write(timg(:,:,ii));
    end
end

t.close();
