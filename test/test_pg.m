
    opt             = struct;
    opt.directory   = '/Users/balbasty/Desktop/model/output';
    opt.fnames.a    = fullfile('/Users/balbasty/Desktop/model', 'IXIcrap_mu.nii');
    opt.fnames.w    = fullfile('/Users/balbasty/Desktop/model', 'IXIcrap_Wv.nii');
    opt.fnames.f    = fullfile('/Users/balbasty/Desktop/model/input', 'IXI002-Guys-0828-T1.img');
    opt.model       = struct('name', 'categorical');
    opt.prm         = s.v_settings;
    opt.regz        = WWv;
    opt.debug       = false;
    opt.gnit        = 20;
    opt.par         = false;

    opt = pg_registration(opt);
    