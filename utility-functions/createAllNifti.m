function createAllNifti(dat, model, opt)

    if isfield(opt, 'prm')
        prm = opt.prm;
    elseif isfield(opt, 'pg') && isfield(opt.pg, 'prm')
        prm = opt.pg.prm;
    else
        prm = [0 0 0 0 0];
    end
    
    if isfield(model, 'Mmu')
        Mmu = model.Mmu;
    elseif isfield(model, 'M')
        Mmu = model.M;
    else
        Mmu = eye(4);
    end
                
    createDat(dat, prm);
    createModel(model, prm, Mmu);

end

function createDat(dat, prm, Mf)

    for n=1:numel(dat)
        if nargin < 3
            if isfield(dat(n), 'Mf')
                Mf = dat(n).Mf;
            elseif isfield(dat(n), 'f') && isfield(dat(n).f, 'M')
                Mf = dat(n).f.M;
            else
                Mf = eye(4);
            end
        end
        fields = fieldnames(dat(n));
        for i=1:numel(fields)
            field = fields{i};
            if isstruct(dat(n).(field))
                createDat(dat(n).(field), prm, Mf);
            else
                if ~strcmpi(field, 'f')
                    switch field
                        case 'v',    descrip = sprintf('Initial velocity (%f %f %f %f %f)', prm);
                        case 'r',    descrip = sprintf('Residual velocity (%f %f %f %f %f)', prm);
                        case 'iphi', descrip = 'Inverse diffeomorphism';
                        case 'phi',  descrip = 'Direct diffeomophism';
                        case 'jac',  descrip = 'Direct Jacobian';
                        case 'ipsi', descrip = 'Inverse transform';
                        case 'psi',  descrip = 'Direct transform';
                        case 'wmu',  descrip = 'Warped template';
                        case 'pf',   descrip = 'Pushed image';
                        case 'c',    descrip = 'Count image';
                        otherwise,   descrip = 'Shape toolbox internal';
                    end
                    if isa(dat(n).(field), 'file_array') ...
                            && checkarray(dat(n).(field)) ...
                            && ~strcmpi(dat(n).(field).permission, 'ro')
                        nii         = nifti;
                        nii.mat0    = Mf;
                        nii.mat     = Mf;
                        nii.descrip = descrip;
                        nii.dat     = dat(n).(field);
                        if numel(nii.dat.dim) > 3
                            nii.dat.dim = [nii.dat.dim(1:3) 1 nii.dat.dim(4:end)];
                        end
                        create(nii);
                    end
                end
            end
        end
    end
    
end

function createModel(model, prm, Mmu)

    fields = fieldnames(model);
    for i=1:numel(fields) 
        field = fields{i};
        if isstruct(model.(field))
            createModel(model.(field), prm, Mmu);
        else
            switch field
                case 'w',    descrip = sprintf('Principal geodesic (%f %f %f %f %f)', prm);
                case 'a',    descrip = 'Log-template';
                case 'mu',   descrip = 'Template';
                case 'gmu',  descrip = 'Template gradients';
                case 'dw',   descrip = 'PG search direction';
                case 'gw',   descrip = 'PG gradient';
                case 'hw',   descrip = 'PG hessian';
                otherwise,   descrip = 'Shape toolbox internal';
            end
            if isa(model.(field), 'file_array') ...
                    && checkarray(model.(field)) ...
                    && ~strcmpi(model.(field).permission, 'ro')
                nii         = nifti;
                nii.mat0    = Mmu;
                nii.mat     = Mmu;
                nii.descrip = descrip;
                nii.dat     = model.(field);
                if numel(nii.dat.dim) > 3
                    nii.dat.dim = [nii.dat.dim(1:3) 1 nii.dat.dim(4:end)];
                end
                create(nii);
            end
        end
    end
    
end