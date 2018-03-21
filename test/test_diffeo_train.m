function test_diffeo_train

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
    opt.dir.model    = '/Users/balbasty/Desktop/IXI_2D/output_diffeo';
    opt.dir.dat      = opt.dir.model;
    opt.model        = struct('name', 'categorical');

    % Special stuff for Holly
    opt.split.par    = 0;
    opt.split.batch  = N;
    
    diffeo_train(input, opt);
    