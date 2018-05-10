function varargout = ghMatchingVel(model, mu, f, varargin)
%__________________________________________________________________________
%
% Compute gradient/hessian of the **negative** log-likelihood of the 
% matching term with respect to changes in the complete velocity.
%
% -------------------------------------------------------------------------
%
% FORMAT [(g), (h, htype)] = ghMatchingVel(model, mu, f, (ga), ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Template image ([mx my mz k])
% f     - Observed image pushed into template space ([mx my mz k])
%
% OPTIONAL
% --------
% ga    - Spatial gradients of the log-probability template.
%         > If ga provided, but not ipsi, compute pointwise gM * gA
%           In this case, ga size is [mx my mz k 3]
%
% KEYWORD ARGUMENTS
% -----------------
% ipsi    - Inverse (subj to template) warp [mx my mz 3]
%           > If provided on top of ga, compute push(gM) * gA
%             In this case, ga size is [nx ny nz k 3]
% lat     - Output lattice size (not needed if ga provided)
% circ    - (Push only) Boundary conditions for pushing [1]
% count   - Pushed voxel count ([nx ny nz nc])
% bb      - Bounding box (if different between template and pushed image)
% hessian - Return only [h, htype] (not g)
% loop    - Specify how to split data processing
%           ('slice', 'component' or 'none' [default])
% par     - If true, parallelise processing. [default: false]
%
% OUTPUT
% ------
% g     - Gradient
%         > If ga and ipsi: [nx ny nz 3]
%         > Else if ga:     [mx my mz 3]
%         > Else if ipsi:   [nx ny nz k]
%         > Else:           [mx my mz k]
% h     - Hessian
%         > If ga and ipsi: [nx ny nz 6]
%         > Else if ga:     [mx my mz 6]
%         > Else if ipsi:   [nx ny nz k(k+1)/2]  for Categorical model
%                           [nx ny nz k]         for other models
%         > Else:           [mx my mz k(k+1)/2]  for Categorical model
%                           [mx my mz k]         for other models
% htype - Shape of the hessian
%         > If ga:          'symtensor'
%         > Else:           'symtensor'   for Categorical model
%                           'diagonal'    for other models
%
%--------------------------------------------------------------------------
%
% Models:
% normal      - Gaussian noise model     (F ~ Mu(IPsi) + N(0,s))
% laplace     - Laplace noise model      (F ~ Mu(IPsi) + L(0,b))
% bernoulli   - True/False realisation   (F ~ Ber(Mu(IPsi)))
% categorical - Multiclass realisation   (F ~ Cat(Mu(IPsi)))
%
% Gradients actually take the form of vector fields ([nx ny nz 3 k]) 
% and  Hessians take the form of tensor fields ([nx ny nz 3 k k]), 
% but the  "vector" dimension (which consists of a pointwise multiplication 
% with  the template spatial gradients) can be performed outside for  
% computational reasons.
%__________________________________________________________________________

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghMatchingVel';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('gmu', []);
    p.addParameter('ipsi',   []);
    p.addParameter('lat',    []);
    p.addParameter('circ',    1);
    p.addParameter('count',  [],     @(X) isnumeric(X) || isa(X, 'file_array'));
    p.addParameter('bb',     struct, @isstruct);
    p.addParameter('hessian', false, @islogical);
    p.addParameter('loop',   '',     @ischar);
    p.addParameter('par',    false,  @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false,  @isscalar);
    p.parse(model, mu, f, varargin{:});
    gmu     = p.Results.gmu;
    ipsi    = p.Results.ipsi;
    lat     = p.Results.lat;
    circ    = p.Results.circ;
    count   = p.Results.count;
    bb      = p.Results.bb;
    hessian = p.Results.hessian;
    loop    = p.Results.loop;
    par     = p.Results.par;
    output  = p.Results.output;
    debug   = p.Results.debug;
    
    if debug, fprintf('* ghMatchingVel\n'); end
    
    
    % --- Compute gradient and hessian (select case)
    switch lower(model.name)
        % E = -log p(f | mu, v)
        case {'normal', 'gaussian', 'l2'}
            if ~isfield(model, 'sigma2')
                model.sigma2 = 1;
            end
            [varargout{1:nargout}] = ...
                ghNormal(mu, f, model.sigma2, gmu, ...
                        'ipsi',    ipsi, ...
                        'lat',     lat, ...
                        'circ',    circ, ...
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
                ghLaplace(mu, f, model.b, gmu, ...
                        'ipsi',    ipsi, ...
                        'lat',     lat, ...
                        'circ',    circ, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'bernoulli', 'binomial', 'binary'}
            [varargout{1:nargout}] = ...
                ghBernoulli(mu, f, gmu, ...
                        'ipsi',    ipsi, ...
                        'lat',     lat, ...
                        'circ',    circ, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'categorical', 'multinomial'}
            [varargout{1:nargout}] = ...
                ghCategorical(mu, f, gmu, ...
                        'ipsi',    ipsi, ...
                        'lat',     lat, ...
                        'circ',    circ, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);
            
        otherwise
            error('Unknown likelihood function.');
    end
    
    
end