function a = reconstructLogProbaTemplate(na, varargin)
% FORMAT a = reconstructLogProbaTemplate(na, ...)
% 
% ** Required **
% na - Encoding of the log-probability template in the null space
% ** Keyword arguments **
% loop - How to split data processing ('none', 'slice') [auto]
% par  - If true, parallelise processing [false]
% ** Output **
% a  - Log-probabilities in the "classical" space
%
% This function allows to deal with log-probabilities which are encoded in
% a subspace of dimension K-1, enforcing the condition that the "actual" 
% K log-probabilities sum to zero. This allows to make the log-probability
% encoding bijective instead of surjective (because probabilities sum to 1,
% they actually belong to a subspace of dimension K-1).

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructLogProbaTemplate';
    p.addRequired('na',  @checkarray);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none'})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(na, varargin{:});
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* reconstructLogProbaTemplate\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(na, 'file_array'), size(na, 3));
    
    % --- Read dim info
    dim = [size(na) 1 1];
    lat = dim(1:3);
    nc  = dim(4)+1;
    
    R = null(ones(1, nc));
    
    % --- Allocate output
    a = prepareOnDisk(output, [lat nc]);
    
    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end;
        a(:,:,:,:) = reshape(reshape(na, [], nc-1) * R', [lat nc]);
        
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
            a(:,:,z,:) = reshape(reshape(na(:,:,z,:), [], nc-1) * R', [latz nc])
        end
        
    end
    
    % --- Write on disk
    a = saveOnDisk(output, a);

end