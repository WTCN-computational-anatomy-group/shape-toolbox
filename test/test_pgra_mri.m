
    fprintf('[Build input file arrays]\n')
    N = 15;
    input_dir = '/Users/balbasty/Devel/ucl/data-balbasty/IXI-2D';
    subjdirs  = spm_select('List', input_dir, 'dir');
    falist = cell(1, N);
    for n=1:N
        subjdir = fullfile(input_dir, deblank(subjdirs(n,:)));
        files   = spm_select('List', subjdir, '\.nii$');
        modalities = cell(1, size(files, 1));
        for i=1:size(files, 1)
            nii = nifti(fullfile(subjdir, deblank(files(i,:))));
            modalities{i} = nii.dat;
        end
        facat = cat(4, modalities{:});
        falist{n} = facat;
    end
    dat = struct('f', falist);
    
    for n=1:N
        f = numeric(dat(n).f);
        msk = ~(f > 0);
        f(msk) = inf;
        dat(n).f(:,:,:,:) = f(:,:,:,:);
    end

    opt              = struct;
    opt.directory    = '/Users/balbasty/Desktop/model/output2';
    opt.model        = struct('name', 'normal', 'sigma2', 1000);
    opt.K            = 10;
    opt.prm          = 100*[0 0.001 0.02 0.0025 0.005];
    opt.emit         = 1000;
    opt.gnit         = 1;
    opt.par          = inf;
    opt.loop         = 'subject';
    opt.debug        = false;
    opt.batch        = 10;
    opt.Mf           = [-1.5 0   0   0;
                         0   1.5 0   0;
                         0   0   1.5 0;
                         0   0   0   1];
    opt.nz0          = 0;
    opt.wpz0         = [1 1];
    opt.affine_basis = affine_basis(12, '2d');
    opt.armijo       = 16;
    opt.gnit         = 5;
    opt.fwhm         = 6;

    fprintf('[Launch pgra_model]\n')
    [model, dat] = pgra_model(opt, dat);
    