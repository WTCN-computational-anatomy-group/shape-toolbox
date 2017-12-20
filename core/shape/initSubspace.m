function [w, ww] = initSubspace(lat, K, varargin)
% FORMAT [w, (ww)] = initSubspace(lat, K,
%                                 ('type', type),
%                                 ('output', output),
%                                 ('debug', debug))
% lat    - Lattice dimension (nx*ny*nz)
% K      - Number of principal components
% type   - Initialisation type ['zero']
% output - array or file_array where to store the result [unused]
% debug  - Print debug info [false]
%
% Init the principal subspace.

    % --- Parser
    p = inputParser;
    p.addRequired('lat');
    p.addRequired('K');
    p.addParameter('type', 'zero', @ischar);
    p.addParameter('output', {[], []});
    p.addParameter('debug', false, @isscalar);
    p.parse(lat, K, varargin{:});
    type   = p.Results.type;
    output = p.Results.output;
    debug  = p.Results.debug;
   
    if debug, fprintf('* initSubspace\n'); end;
     
    % --- Check output size
    if ~iscell(output)
        output = {output};
    end
    for i=numel(output)+1:nargout
        output{i} = [];
    end
    
    % --- Prepare output
    w = prepareOnDisk(output{1}, [lat 3 K]);
    if nargout > 1
        ww = prepareOnDisk(output{2}, [K K]);
    end
    
    % --- Initialise
    switch lower(type)
        case 'zero'
            for k=1:size(w, 5)
                w(:,:,:,:,k) = 0;
            end
            if nargout > 2
                ww(:,:) = eye(K);
            end
        otherwise
            error('Unknown initialisation type ''%s''', type);
    end

end