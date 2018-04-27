function na = rotateLogTemplate(a, varargin)
% FORMAT na = rotateLogTemplate(na, ...)
% 
% REQUIRED
% --------
% a - Log-probability template in the "classical" space
%
% KEYWORD ARGUMENTS
% -----------------
% loop   - How to split data processing ('none', 'slice') [auto]
% par    - If true, parallelise processing [false]
% output - File array to store the output [not used]
% debug  - Debugging talk [false]
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
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none'})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(a, varargin{:});
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* rotateLogTemplate\n'); end
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(a, 'file_array'), size(a, 3));
    
    % --- Read dim info
    dim = [size(a) 1 1];
    lat = dim(1:3);
    nc  = dim(4);
    
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
    na = prepareOnDisk(output, [lat nc-1]);
    
    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end
        na(:,:,:,:) = reshape(reshape(a, [], nc) * R, [lat nc-1]);
        
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
        latz = [lat(1:2) 1];
        parfor (z=1:lat(3), par)
            na(:,:,z,:) = reshape(reshape(a(:,:,z,:), [], nc) * R, [latz nc-1])
        end
        
    end
    
    % --- Write on disk
    na = saveOnDisk(output, na);
    if tmp_file
        delete(fname);
    end

end