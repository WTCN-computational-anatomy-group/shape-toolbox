function [opt, dat] = pg_model_input(opt, dat)
% FORMAT [opt, dat] = pg_model_input(opt, dat)
%
% Parse provided observed images and start setting the data structure.
% Should be called first.

    % ---------------------------------------------------------------------
    % Input images
    % ---------------------------------------------------------------------
    
    if ~isfield(opt, 'fnames')
        opt.fnames = struct;
        if ~isfield(opt.fnames, 'dat')
            opt.fnames.dat = struct;
        end
    end
    
    % --- Filenames provided
    if ~isfield(dat, 'f')
        
        % - If needed convert from char array to cell array
        if isfield(opt.fnames.dat, 'f') ...
                && ischar(opt.fnames.dat.f) ...
                && size(opt.fnames.dat.f, 1) > 1
            fnames = cell(1, size(opt.fnames.dat.f, 1));
            for n=1:numel(fnames)
                fnames{n} = deblank(opt.fnames.dat.f(n,:));
            end
            opt.fnames.dat.f = fnames;
        end
        % - Check there is an input
        if ~isfield(opt.fnames.dat, 'f') || ~iscell(opt.fnames.dat.f)
            error('A set of observed images must be provided')
        end
        % - Build data structure
        f  = cell(1, numel(opt.fnames.dat.f));
        Mf = cell(1, numel(opt.fnames.dat.f));
        for i=1:numel(opt.fnames.dat.f)
            if ~exist(opt.fnames.dat.f{i}, 'file')
                warning('Image %s does not exist.', opt.fnames.dat.f{i})
                f  = f(1:end-1);
                Mf = Mf(1:end-1);
            else
                n = nifti(opt.fnames.dat.f{i});
                f{i}  = n.dat;
                Mf{i} = n.mat0;
            end
        end
        dat = struct('f', f, 'Mf', Mf);
        clear f Mf
        if isfield(opt, 'Mf')
            [dat.Mf] = deal(opt.Mf);
        end
        
    % --- Data provided (file_array or array)
    else
        
        if isfield(opt, 'Mf')
            [dat.Mf] = deal(opt.Mf);
            
        else
            for i=1:numel(dat)
                if isa(dat(i).f, 'file_array')
                    fnames = {dat(i).f.fname};
                    if endsWith(fnames{1}, '.nii')
                        n = nifti(fnames{1});
                        dat(i).Mf = n.mat0;
                    else
                        dat(i).Mf = eye(4);
                    end
                else
                    dat(i).Mf = eye(4);
                end
            end
            
        end
        
    end
    opt.N = numel(dat);

end