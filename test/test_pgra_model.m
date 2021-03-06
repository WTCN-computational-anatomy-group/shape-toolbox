
    N = 100;
    input_dir = '/Users/balbasty/Desktop/IXI_2D/input';
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

    opt = struct;
    
    % General options
    opt.dir.model    = '/Users/balbasty/Desktop/IXI_2D/output_pgra';
    opt.dir.dat      = opt.dir.model;
    opt.model        = struct('name', 'categorical');
    opt.pg.K         = 19;  % Change / number of subjects

    % Special stuff for Holly
    opt.split.par    = 0;
    opt.split.batch  = N;
    opt.ondisk.dat.v.g = true;
    opt.ondisk.dat.v.h = true;
    opt.ui.ftrack = true;
    
%     opt              = struct;
%     opt.dir.model    = '/Users/balbasty/Desktop/model/output';
%     opt.dir.dat      = opt.dir.model;
%     opt.model        = struct('name', 'categorical');
%     opt.pg.prm       = [0 0.001 0.02 0.0025 0.005];
%     opt.split.par    = 0;
%     opt.pg.K         = 19;
%     opt.lb.moving    = 1;
%     opt.optimise.q   = false;

    [model, dat] = pgra_model(input, opt);
    