function createAllNifti(dat, model, opt)

    % ---
    % DAT
    % ---
    fields = fieldnames(dat);
    for i=1:numel(fields)
        field = fields{i};
        if ~strcmpi(field, 'f')
            switch field
                case 'v',    descrip = sprintf('Initial velocity (%f %f %f %f %f)', opt.prm);
                case 'r',    descrip = sprintf('Residual velocity (%f %f %f %f %f)', model.lambda*opt.prm);
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
            for n=1:numel(dat)
                if isa(dat(n).(field), 'file_array') && checkarray(dat(n).(field))
                    nii         = nifti;
                    nii.mat0    = dat(n).Mf;
                    nii.mat     = dat(n).Mf;
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

    % -----
    % MODEL
    % -----
    fields = fieldnames(model);
    for i=1:numel(fields)
        field = fields{i};
        switch field
            case 'w',    descrip = sprintf('Principal geodesic (%f %f %f %f %f)', opt.prm);
            case 'a',    descrip = 'Log-template';
            case 'mu',   descrip = 'Template';
            case 'gmu',  descrip = 'Template gradients';
            case 'dw',   descrip = 'PG search direction';
            case 'gw',   descrip = 'PG gradient';
            case 'hw',   descrip = 'PG hessian';
            otherwise,   descrip = 'Shape toolbox internal';
        end
        if isa(model.(field), 'file_array') && checkarray(model.(field))
            nii         = nifti;
            nii.mat0    = model.Mmu;
            nii.mat     = model.Mmu;
            nii.descrip = descrip;
            nii.dat     = model.(field);
            if numel(nii.dat.dim) > 3
                nii.dat.dim = [nii.dat.dim(1:3) 1 nii.dat.dim(4:end)];
            end
            create(nii);
        end
    end
end