function varargout = ghTemplate(model, mu, f, varargin)
%__________________________________________________________________________
%
% Compute gradient/hessian of the **negative** log-likelihood of the 
% matching term with respect to changes in the template values.
%
% -------------------------------------------------------------------------
%
% FORMAT [(g), (h, htype)] = ghTemplate(model, mu, f, ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Template image ([mx my mz k]/[nx ny nz k])
% f     - Observed image ([mx my mz k]/[nx ny nz k])
%
% KEYWORD ARGUMENTS
% -----------------
% ipsi    - Inverse (subj to template) warp [mx my mz 3]
%           > If provided, compute push(gM)
% circ    - (Push only) Boundary conditions for pushing [1]
% lat     - Output lattice size ([mx my mz])
% count   - Pushed voxel count ([nx ny nz])
% bb      - Bounding box (if different between template and pushed image)
% hessian - Return only [h, htype] (not g)
% rotate  - [categorical only] rotate out null space [false]
% loop    - Specify how to split data processing
%           ('slice', 'component' or 'none' [default])
% par     - If true, parallelise processing. [default: false]
%
% OUTPUT
% ------
% g     - Gradient
%         > If ipsi:        [nx ny nz k]
%         > Else:           [mx my mz k]
% h     - Hessian
%         > If ipsi:        [nx ny nz k(k+1)/2]  for Categorical model
%                           [nx ny nz k]         for other models
%         > Else:           [mx my mz k(k+1)/2]  for Categorical model
%                           [mx my mz k]         for other models
% htype - Shape of the hessian
%         > 'symtensor'   for Categorical model
%           'diagonal'    for other models
%
%--------------------------------------------------------------------------
%
% Models:
% normal      - Gaussian noise model     (F ~ Mu(IPsi) + N(0,s))
% laplace     - Laplace noise model      (F ~ Mu(IPsi) + L(0,b))
% bernoulli   - True/False realisation   (F ~ Ber(Mu(IPsi)))
% categorical - Multiclass realisation   (F ~ Cat(Mu(IPsi)))
%
% Gradients actually take the form of vector fields ([nx ny nz k]) 
% and  Hessians take the form of tensor fields ([nx ny nz k k]).
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghTemplate';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addParameter('ipsi',   []);
    p.addParameter('circ',     1);
    p.addParameter('lat',    []);
    p.addParameter('count',  [],     @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('rotate',  false, @islogical);
    p.addParameter('hessian', false, @islogical);
    p.addParameter('loop',   '',     @ischar);
    p.addParameter('par',    false,  @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false,  @isscalar);
    p.parse(model, mu, f, varargin{:});
    ipsi    = p.Results.ipsi;
    circ    = p.Results.circ;
    lat     = p.Results.lat;
    count   = p.Results.count;
    bb      = p.Results.bb;
    rotate  = p.Results.rotate;
    hessian = p.Results.hessian;
    loop    = p.Results.loop;
    par     = p.Results.par;
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghTemplate\n'); end
    
    
    % --- Compute gradient and hessian (select case)
    switch lower(model.name)
        % E = -log p(f | mu, v)
        case {'normal', 'gaussian', 'l2'}
            if ~isfield(model, 'sigma2')
                model.sigma2 = 1;
            end
            [varargout{1:nargout}] = ...
                ghNormal(mu, f, model.sigma2, ...
                        'ipsi',    ipsi, ...
                        'circ',    circ, ...
                        'lat',     lat, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'laplace', 'l1'}
            if ~isfield(model, 'b')
                model.b = 1;
            end
            [varargout{1:nargout}] = ...
                ghLaplace(mu, f, model.b, ...
                        'ipsi',    ipsi, ...
                        'circ',    circ, ...
                        'lat',     lat, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'bernoulli', 'binomial', 'binary'}
            [varargout{1:nargout}] = ...
                ghBernoulli(mu, f, ...
                        'ipsi',    ipsi, ...
                        'circ',    circ, ...
                        'lat',     lat, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'categorical', 'multinomial'}
            [varargout{1:nargout}] = ...
                ghCategorical(mu, f, ...
                        'ipsi',    ipsi, ...
                        'circ',    circ, ...
                        'lat',     lat, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'rotate',  rotate, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);
            
        otherwise
            error('Unknown likelihood function.');
    end
    
    
end