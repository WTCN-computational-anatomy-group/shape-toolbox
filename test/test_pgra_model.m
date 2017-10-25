
    N = 250;
    input_dir = '/Users/balbasty/Desktop/model/input';
    files     = spm_select('List', input_dir, '\.img$');
    if ~isfinite(N)
        N         = size(files, 1);
    end
    fnames    = cell(1, N);
    for n=1:N
        fnames{n} = fullfile(input_dir, deblank(files(n,:)));
    end
    

    opt              = struct;
    opt.directory    = '/Users/balbasty/Desktop/model/output';
    opt.fnames.dat.f = fnames;
    opt.model        = struct('name', 'categorical');
    opt.K            = 10;
    opt.prm          = [0 0.001 0.02 0.0025 0.005];
    opt.emit         = 1000;
    opt.gnit         = 1;
    opt.par          = inf;
    opt.loop         = 'subject';
    opt.debug        = false;
    opt.batch        = 10;
    opt.nz0          = 0;
    opt.wpz0         = [1 1];
    opt.fwhm         = 6;
    opt.affine_basis = affine_basis(12, '2d');

    [model, dat] = pgra_model(opt);
    