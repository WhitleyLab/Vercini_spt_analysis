function imwritestack(fname,imS);
%function imwritestack(fname,imS);
%assumes tif
fnameFinal = fname;
fname = [tempdir,filesep(),'tempIm.tif'];
nFr = size(imS,3);
for ii = 1:nFr
    if ii == 1
        writeMode = 'overwrite';
    else
        writeMode = 'append';
    end
    imwrite(imS(:,:,ii),fname,'WriteMode',writeMode);
end
movefile(fname,fnameFinal);
