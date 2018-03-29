function [mu, s] = reconstructProbaTemplate(a, varargin)
% _________________________________________________________________________
%
% Reconstruct the template probability maps from their log-space. 
% Only useful for bernoulli and categorical matching terms.
%
% -------------------------------------------------------------------------
% FORMAT mu = reconstructProbaTemplate(a, ...)
%
% REQUIRED
% --------
% a    - Log-Template
%
% KEYWORD ARGUMENTS
% -----------------
% type   - Mapping type ('sigmoid' or 'softmax')                    [auto]
% scale  - Image-wide probability scaling factor                    [ones]
% loop   - Specify how to split data processing ('slice'/'none')    [auto]
% par    - If true, parallelise processing                          [false]
% output - file_array where to store the output                     []
% debug  - Debugging talk                                           [false]
%
% OUTPUT
% ------
% mu  - Reconstructed template
% s   - Sufficient statistics for scale update, i.e.:
%         s_k = sum[ exp(a_k)/sum[exp(a_l + scl_l)] ]
% _________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructProbaTemplate';
    p.addRequired('a',              @checkarray);
    p.addParameter('type',   '',    @(X) ischar(X) && any(strcmpi(X, {'sigmoid', 'softmax'})));
    p.addParameter('scale',  1,     @isnumeric);
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(a, varargin{:});
    type   = p.Results.type;
    scale  = p.Results.scale;
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* reconstructProbaTemplate\n'); end

    % --- Automatic par/loop setting
    [par, loop] = autoParLoop(par, loop, isa(a, 'file_array'), size(a, 3), 1);
    
    % --- Read dimensions
    dim = [size(a) 1 1];
    dim = dim(1:4);
    lat = dim(1:3);  % Template lattice dimension
    nc  = dim(4);    % Number of classes/modalities
    
    if isempty(type)
        if nc == 1
            type = 'sigmoid';
        else
            type = 'softmax';
        end
    end
    
    % --- Scale factor
    if numel(scale) == 1
        scale = repmat(scale, 1, nc);
    end
    scale = reshape(scale, [1 1 1 nc]);
  
    % --- Reserve space on disk
    mu = prepareOnDisk(output, [lat nc]);
    s  = 0;
    if strcmpi(type, 'sigmoid')
        s = nan;
    end

    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end
        if strcmpi(type, 'sigmoid')
            mu(:,:,:,:) = sigmoid(a);
        else
            mu(:,:,:,:) = softmax(a, scale);
            s = s + sum(sum(sum(softmax(a),1),2),3);
        end
        
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
        if strcmpi(type, 'sigmoid')
            parfor (z=1:lat(3), par)
                mu(:,:,z,:) = sigmoid(a(:,:,z,:));
            end
        else
            parfor (z=1:lat(3), par)
                mu(:,:,z,:) = softmax(a(:,:,z,:), scale);
                s = s + sum(sum(sum(softmax(a(:,:,z,:)),1),2),3);
            end
        end
        
    end
    
    s = s(:);

    % --- Write on disk
    mu = saveOnDisk(output, mu);
        
end

function mu = sigmoid(a)
    a = exp(single(numeric(a)));
    mu = a./(1 + a);
end

function a = softmax(a, scale)
    if nargin < 2
        scale = 1;
    end
    a = single(numeric(a));               % read from disk if needed
    a = bsxfun(@plus, a, log(scale));     % rescale probabilities
    a = bsxfun(@minus, a, max(a, [], 4)); % safe softmax -> avoid overflow
    a = exp(a);                           % exponentiate
    a = bsxfun(@rdivide, a, sum(a, 4));   % normalise
end