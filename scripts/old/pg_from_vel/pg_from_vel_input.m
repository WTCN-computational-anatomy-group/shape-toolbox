function [opt, dat] = pg_from_vel_input(opt, dat)
% FORMAT [opt, dat] = pg_from_vel_input(opt, dat)
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
    if ~isfield(dat, 'v')
        
        % - If needed convert from char array to cell array
        if isfield(opt.fnames.dat, 'v') ...
                && ischar(opt.fnames.dat.v) ...
                && size(opt.fnames.dat.v, 1) > 1
            fnames = cell(1, size(opt.fnames.dat.v, 1));
            for n=1:numel(fnames)
                fnames{n} = deblank(opt.fnames.dat.v(n,:));
            end
            opt.fnames.dat.v = fnames;
        end
        % - Check there is an input
        if ~isfield(opt.fnames.dat, 'v') || ~iscell(opt.fnames.dat.v)
            error('A set of observed images must be provided')
        end
        % - Build data structure
        v  = cell(1, numel(opt.fnames.dat.v));
        Mf = cell(1, numel(opt.fnames.dat.v));
        for i=1:numel(opt.fnames.dat.v)
            if ~exist(opt.fnames.dat.v{i}, 'file')
                warning('Image %s does not exist.', opt.fnames.dat.v{i})
                v  = v(1:end-1);
                Mf = Mf(1:end-1);
            else
                n = nifti(opt.fnames.dat.v{i});
                v{i}  = n.dat;
                Mf{i} = n.mat0;
            end
        end
        dat = struct('v', v, 'Mf', Mf);
        clear v Mf
        if isfield(opt, 'Mf')
            [dat.Mf] = deal(opt.Mf);
        end
        
    % --- Data provided (file_array or array)
    else
        
        if isfield(opt, 'Mf')
            [dat.Mf] = deal(opt.Mf);
            
        else
            for i=1:numel(dat)
                if isa(dat(i).v, 'file_array')
                    fnames = {dat(i).v.fname};
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