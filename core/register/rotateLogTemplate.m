function na = rotateLogTemplate(a, varargin)
% FORMAT na = rotateLogTemplate(na, ...)
% 
% REQUIRED
% --------
% a - Log-probability template in the "classical" space
%
% KEYWORD ARGUMENTS
% -----------------
% hessian - The input is a Hessian and must be rotated on both sides [false]
% loop    - How to split data processing ('none', 'slice') [auto]
% par     - If true, parallelise processing [false]
% output  - File array to store the output [not used]
% debug   - Debugging talk [false]
%
% OUTPUT
% ------
% na  - Log-probabilities in the null space
%
% This transform projects log-probabilities to a subspace of dimension K-1
% enforcing the condition that the "actual" K log-probabilities sum to zero.
% This allows to make the log-probability encoding bijective instead of 
% surjective (because probabilities sum to 1, they actually belong to a 
% subspace of dimension K-1) and makes the Hessian more stable.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'rotateLogTemplate';
    p.addRequired('a',  @checkarray);
    p.addParameter('hessian', false,@(X) isscalar(X) && islogical(X));
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none'})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(a, varargin{:});
    hessian = p.Results.hessian;
    par     = p.Results.par;
    loop    = p.Results.loop;
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* rotateLogTemplate\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(a, 'file_array'), size(a, 3));
    
    % --- Read dim info
    dim = [size(a) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
    % --- Hessian case
    if hessian
        sparse = numel(size(a)) < 5;
        if sparse
            K = nc;
            [~, nc] = spm_matcomp('SymIndices', K);
        end
    end
    
    R = null(ones(1, nc));
    
    % --- Check input == output
    tmp_file = false;
    if isa(a, 'file_array') && isa(output, 'file_array') ...
            && strcmpi(a.fname, output.fname)
        [path, name, ext] = fileparts(a.fname);
        fname = fullfile(path, ['tmp_' name ext]);
        copyfile(a.fname, fname);
        a.fname  = fname;
        tmp_file = true;
    end
    
    % --- Allocate output
    if hessian
        if sparse
            na = prepareOnDisk(output, [lat nc*(nc-1)/2]);
        else
            na = prepareOnDisk(output, [lat nc-1 nc-1]);
        end
    else
        na = prepareOnDisk(output, [lat nc-1]);
    end
    
    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end
        na(:,:,:,:,:) = rotate(R, a, hessian);
        
    % --- Loop on slices
    elseif strcmpi(loop, 'slice')
        % - Parallelised
        if debug
            if par > 0
                fprintf('   - Parallelise on slices\n'); 
            else
                fprintf('   - Serialise on slices\n'); 
            end
        end
        if par > 0
            if isa(a, 'file_array')
                parfor (z=1:lat(3), par)
                    na(:,:,z,:,:) = rotate(R, slicevol(a, z, 3), hessian);
                end
            else
                parfor (z=1:lat(3), par)
                    na(:,:,z,:,:) = rotate(R, a(:,:,z,:,:), hessian);
                end
            end
        else
            for z=1:lat(3)
                na(:,:,z,:,:) = rotate(R, a(:,:,z,:,:), hessian);
            end
        end
        
    end
    
    % --- Write on disk
    na = saveOnDisk(output, na);
    if tmp_file
        delete(fname);
    end

end

function a = rotate(R, a, hessian)
    dim = [size(a) 1 1];
    lat = dim(1:3);
    N   = dim(4);
    a   = numeric(a);
    if hessian 
        sparse = numel(size(a)) < 5;
        if sparse
            K = size(a, 4);
            [ind, N] = spm_matcomp('SymIndices', K);
            tmp = a;
            a = zeros([lat N N], 'single');
            for i=1:N
                for j=1:N
                    a(:,:,:,i,j) = tmp(:,:,:,ind(i,j));
                end
            end
            clear tmp
        end
        a = reshape(reshape(a, [], N) * R, [lat N N-1]);
        a = shiftdim(a, 3);
        a = reshape(R' * reshape(a, N, []), [N-1 N-1 lat]);
        a = shiftdim(a, 2);
        a = reshape(a, [lat N-1 N-1]); % in case there were singleton dimensions
        if sparse
            tmp = a;
            [ind, K] = spm_matcomp('SymIndices', N-1, 'n');
            a = zeros([lat K], 'single');
            for i=1:N-1
                a(:,:,:,ind(i,i)) = tmp(:,:,:,i,i);
                for j=i+1:N-1
                    a(:,:,:,ind(i,j)) = tmp(:,:,:,i,j);
                end
            end
            clear tmp
        end
    else
        a = reshape(reshape(a, [], N) * R, [lat N-1]);
    end
end