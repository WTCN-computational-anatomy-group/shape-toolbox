function mu = reconstructProbaTemplate(a, varargin)
% FORMAT mu = reconstructProbaTemplate(a, ('type', type),
%                                      ('loop', loop), ('par', par))
%
% ** Required **
% a    - Log-Template
% ** Keyword arguments **
% type - Mapping type ('sigmoid' or 'softmax') [auto]
% loop - Specify how to split data processing ('slice' or 'none') [auto]
% par  - If true, parallelise processing [false]
% ** Output **
% mu  - Reconstructed template
%
% Reconstruct the template probability maps from their log-space. Only
% useful for bernoulli and categorical matching terms.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'reconstructProbaTemplate';
    p.addRequired('a',              @checkarray);
    p.addParameter('type',   '',    @(X) ischar(X) && any(strcmpi(X, {'sigmoid', 'softmax'})));
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(a, varargin{:});
    type   = p.Results.type;
    par    = p.Results.par;
    loop   = p.Results.loop;
    output = p.Results.output;
    debug  = p.Results.debug;
    
    if debug, fprintf('* reconstructProbaTemplate\n'); end;

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
  
    % --- Reserve space on disk
    mu = prepareOnDisk(output, [lat nc]);

    % --- No loop
    if strcmpi(loop, 'none')
        if debug, fprintf('   - No loop\n'); end;
        if strcmpi(type, 'sigmoid')
            mu(:,:,:,:) = sigmoid(a);
        else
            mu(:,:,:,:) = softmax(a);
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
                mu(:,:,z,:) = softmax(a(:,:,z,:));
            end
        end
        
    end

    % --- Write on disk
    mu = saveOnDisk(output, mu);
        
end

function mu = sigmoid(a)
    a = exp(single(numeric(a)));
    mu = a./ (1 + a);
end

function mu = softmax(a)
    a = exp(single(numeric(a)));
    s = sum(a, 4);
    mu = bsxfun(@rdivide, a, s);
end