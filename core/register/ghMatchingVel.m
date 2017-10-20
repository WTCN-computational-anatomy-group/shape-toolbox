function varargout = ghMatchingVel(model, mu, f, c, varargin)
% FORMAT [g,h] = ghMatchingVel(model, mu, f, c, (gmu), ...)
%
% ** Required **
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu    - Template image ([nx ny nz nc])
% f     - Observed image pushed into template space ([nx ny nz nc])
% c     - Pushed voxel count ([nx ny nz nc])
% ** Optional **
% gmu   - Template spatial gradients.
%         If not provided, do not use: it's deal with outside (useful
%         when computing grad/hess w.r.t. Z)
% ** Keyword arguments **
% loop  - Specify how to split data processing
%         ('slice', 'component' or 'none' [default])
% par   - If true, parallelise processing. [default: false]
% ** Output **
% g     - First derivatives w.r.t. full velocity ([nx ny nz nc])
% h     - Second derivatives w.r.t. full velocity
%         (diagonal approximation: [nx ny nz nc], except for the multinomial 
%         case where it is a symmetric tensor field of virtual size 
%         [nx ny nz nc nc])
% htype - Type of hessian approximation: 'diagonal' or 'symtensor'
%
% Models:
% normal      - Gaussian noise model     (F ~ Mu(IPhi) + N(0,s))
% laplace     - Laplace noise model      (F ~ Mu(IPhi) + L(0,b))
% bernoulli   - True/False realisation   (F ~ Ber(Mu(IPhi)))
% categorical - Multiclass realisation   (F ~ Cat(Mu(IPhi)))
%
% Compute gradient/hessian of the **negative** log-likelihood of the 
% matching term with respect to changes in the complete velocity.
% Gradients actually take the form of vector fields ([nx ny nz ngrad nc]) 
% and  Hessians take the form of tensor fields ([nx ny nz ngrad nc nc]), 
% but the  "vector" dimension (which consists of a pointwise multiplication 
% with  the template spatial gradients) can be performed outside for  
% computational reasons.

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghMatchingVel';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addRequired('c',   @checkarray);
    p.addOptional('gmu', []);
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(model, mu, f, c, varargin{:});
    gmu = p.Results.gmu;
    
    if p.Results.debug, fprintf('* ghMatchingVel\n'); end;
    
    
    % --- Compute gradient and hessian (select case)
    switch lower(model.name)
        % E = -log p(f | mu, v)
        case {'normal', 'gaussian', 'l2'}
            if ~isfield(model, 'sigma2')
                model.sigma2 = 1;
            end
            [varargout{1:nargout}] = ...
                ghNormal(mu, f, c, model.sigma2, gmu, ...
                        'loop',   p.Results.loop, ...
                        'par',    p.Results.par, ...
                        'output', p.Results.output, ...
                        'debug',  p.Results.debug);

        case {'laplace', 'l1'}
            if ~isfield(model, 'b')
                model.b = 1;
            end
            [varargout{1:nargout}] = ...
                ghLaplace(mu, f, c, model.b, gmu, ...
                        'loop',   p.Results.loop, ...
                        'par',    p.Results.par, ...
                        'output', p.Results.output, ...
                        'debug',  p.Results.debug);

        case {'bernoulli', 'binomial', 'binary'}
            [varargout{1:nargout}] = ...
                ghBernoulli(mu, f, c, gmu, ...
                        'loop',   p.Results.loop, ...
                        'par',    p.Results.par, ...
                        'output', p.Results.output, ...
                        'debug',  p.Results.debug);

        case {'categorical', 'multinomial'}
            [varargout{1:nargout}] = ...
                ghCategorical(mu, f, c, gmu, ...
                        'loop',   p.Results.loop, ...
                        'par',    p.Results.par, ...
                        'output', p.Results.output, ...
                        'debug',  p.Results.debug);
            
        otherwise
            error('Unknown likelihood function.');
    end
end