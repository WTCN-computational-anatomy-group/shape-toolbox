function varargout = ghMatchingVel(model, mu, f, varargin)
%__________________________________________________________________________
%
% Compute gradient/hessian of the **negative** log-likelihood of the 
% matching term with respect to changes in the complete velocity.
%
% -------------------------------------------------------------------------
%
% FORMAT [(g),(h),(htype)] = ghMatchingVel(model, mu, f, (ga), (ipsi), ...)
%
% REQUIRED
% --------
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Template image ([nx ny nz nc])
% f     - Observed image pushed into template space ([nx ny nz nc])
%
% OPTIONAL
% --------
% ga    - Spatial gradients of the log-probability template.
%         > If ga provided, but not ipsi, compute pointwise gM * gA
%           In this case, ga size is [mx my mz k]
% ipsi  - Inverse (subj to template) warp [mx my mz 3]
%         > If provided on top of ga, compute push(gM) * gA
%           In this case, ga size is [nx ny nz k]
%
% KEYWORD ARGUMENTS
% -----------------
% count - Pushed voxel count ([nx ny nz nc])
% bb    - Bounding box (if different between template and pushed image)
% hessian - Return only [h, htype] (not g)
% loop    - Specify how to split data processing
%           ('slice', 'component' or 'none' [default])
% par     - If true, parallelise processing. [default: false]
%
% OUTPUT
% ------
% g     - First derivatives w.r.t. full velocity ([nx ny nz nc])
% h     - Second derivatives w.r.t. full velocity
%         (diagonal approximation: [nx ny nz nc], except for the multinomial 
%         case where it is a symmetric tensor field of virtual size 
%         [nx ny nz nc nc])
% htype - Type of hessian approximation: 'diagonal' or 'symtensor'
%
%--------------------------------------------------------------------------
%
% Models:
% normal      - Gaussian noise model     (F ~ Mu(IPhi) + N(0,s))
% laplace     - Laplace noise model      (F ~ Mu(IPhi) + L(0,b))
% bernoulli   - True/False realisation   (F ~ Ber(Mu(IPhi)))
% categorical - Multiclass realisation   (F ~ Cat(Mu(IPhi)))
%
% Gradients actually take the form of vector fields ([nx ny nz ngrad nc]) 
% and  Hessians take the form of tensor fields ([nx ny nz ngrad nc nc]), 
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
    p.addOptional('ipsi', []);
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
                ghNormal(mu, f, model.sigma2, gmu, ipsi, ...
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
                ghLaplace(mu, f, model.b, gmu, ipsi, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'bernoulli', 'binomial', 'binary'}
            [varargout{1:nargout}] = ...
                ghBernoulli(mu, f, gmu, ipsi, ...
                        'count',   count, ...
                        'bb',      bb, ...
                        'hessian', hessian, ...
                        'loop',    loop, ...
                        'par',     par, ...
                        'output',  output, ...
                        'debug',   debug);

        case {'categorical', 'multinomial'}
            [varargout{1:nargout}] = ...
                ghCategorical(mu, f, gmu, ipsi, ...
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