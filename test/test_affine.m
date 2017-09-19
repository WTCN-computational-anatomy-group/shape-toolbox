
    opt             = struct;
    opt.directory   = '/Users/balbasty/Desktop/model/output';
    opt.fnames.a    = fullfile('/Users/balbasty/Desktop/model', 'IXIcrap_mu.nii');
    opt.fnames.f    = fullfile('/Users/balbasty/Desktop/model/input', 'IXI002-Guys-0828-T1.img');
    opt.model       = struct('name', 'categorical');
    opt.debug       = false;
    opt.gnit        = 20;
    opt.affine_basis = affine_basis(12, '2d');
%     opt.regq        = 1E5 * eye(6);
    opt.happrox     = true;
    opt.par = false;
    
    opt = affine_registration(opt);
    