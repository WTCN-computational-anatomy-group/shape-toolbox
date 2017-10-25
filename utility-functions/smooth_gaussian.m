function vol = smooth_gaussian(vol, fwhm)
% FORMAT img = smooth_gaussian(img, fwhm)
% vol  - 2D/3D image to smooth
% fwhm - Full width at half-maximum of the Gaussian kernel
%
% Gaussian smoothing of an image

    if nargin < 2
        fwhm = 0.5;
    end

    if ~fwhm
        return;
    end
    
    if length(fwhm) < 3
        fwhm = padarray(fwhm, 3 - length(fwhm), 'replicate', 'post');
    end
    
    K = size(vol, 4);
    
    lim = ceil(2*fwhm);
    x   = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x = x/sum(x);
    y   = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y = y/sum(y);
    z   = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z = z/sum(z);
    i   = (length(x) - 1)/2;
    j   = (length(y) - 1)/2;
    k   = (length(z) - 1)/2;

    for k1=1:K
        img = vol(:,:,:,k1);
        spm_conv_vol(img,img,x,y,z,-[i j k]);
        vol(:,:,:,k1) = img;
    end
end