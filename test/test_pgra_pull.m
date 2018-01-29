
    N = 20;
    input_dir = '/Users/balbasty/Desktop/model/input';
    files     = spm_select('List', input_dir, '\.img$');
    if ~isfinite(N)
        N         = size(files, 1);
    end
    fnames    = cell(1, N);
    for n=1:N
        fnames{n} = fullfile(input_dir, deblank(files(n,:)));
    end
    input   = struct;
    input.f = fnames;

    opt              = struct;
    opt.dir.model    = '/Users/balbasty/Desktop/model/output';
    opt.dir.dat      = opt.dir.model;
    opt.model        = struct('name', 'categorical');
    opt.pg.prm       = [0 0.001 0.02 0.0025 0.005];
    opt.split.par    = 0;
    opt.pg.K         = 19;
    opt.optimise.pg.w = false;
    opt.optimise.z.z  = false;
    opt.optimise.z.A  = false;
    opt.optimise.q.q  = false;
    opt.optimise.q.A  = false;
    opt.optimise.r.l  = false;

    [model, dat] = pgra_model(input, opt);
    