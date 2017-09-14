function ll = llMatching(model, mu, f, varargin)
% FORMAT ll = llMatching(model, mu, f, (c),
%                        ('type', type), ('loop', loop), ('par', par))
% ** Required **
% model - Structure with fields:
%           * 'name'    : 'normal', 'laplace', 'bernoulli' or 'categorical'
%           * ('sigma2'): Normal variance  [1]
%           * ('b')     : Laplace variance [1]
% mu     - Template image aligned with f
% f      - Observed image aligned with mu
% ** Optional **
% c      - Pushed voxel count.
%            If provided: log-likelihood is computed as in template space
%            Else       : log-likelihood is computed as in image space
% ** Keyword arguments **
% type   - Template type in the bernoulli/categorical case,
%          can be 'proba', 'log' or 'null' [default: 'proba']
% loop   - Specify how to split data processing
%          ('slice', 'component' or 'none' [default: auto])
% par    - If true, parallelise processing. [default: false]
% ** Ouptut **
% ll     - Log-likelihood of the matching term

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'llMatching';
    p.addRequired('model',  @(X) isstruct(X) && isfield(X, 'name'));
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addOptional('c', [], @checkarray);
    p.addParameter('type', 'proba',  @(X) ischar(X) && any(strcmpi(X, {'log', 'proba', 'null'})));
    p.addParameter('loop',   '',    @(X) ischar(X) && any(strcmpi(X, {'slice', 'component', 'none', ''})));
    p.addParameter('par',    false, @isscalar);
    p.addParameter('debug',  false, @isscalar);
    p.parse(model, mu, f, varargin{:});
    c      = p.Results.c;
    type   = p.Results.type;
    loop   = p.Results.loop;
    par    = p.Results.par;
    debug  = p.Results.debug;
    
    if debug, fprintf('* llMatching\n'); end;
    
    % --- Select case
    switch lower(model.name)
        case {'normal', 'gaussian', 'l2'}
            if ~isfield(model, 'sigma2')
                model.sigma2 = 1;
            end
            if isempty(c)
                ll = llNormalWarped(mu, f, model.sigma2, 'loop', loop, 'par', par, 'debug', debug);
            else
                ll = llNormalPushed(mu, f, c, model.sigma2, 'loop', loop, 'par', par, 'debug', debug);
            end
            
        case {'laplace', 'l1'}
            if ~isfield(model, 'b')
                model.b = 1;
            end
            if isempty(c)
                ll = llLaplaceWarped(mu, f, model.b, 'loop', loop, 'par', par, 'debug', debug);
            else
                ll = llLaplacePushed(mu, f, c, model.b, 'loop', loop, 'par', par, 'debug', debug);
            end
            
        case {'bernoulli', 'binomial'}
            if isempty(c)
                ll = llBernoulliWarped(mu, f, type, 'loop', loop, 'par', par, 'debug', debug);
            else
                ll = llBernoulliPushed(mu, f, c, type, 'loop', loop, 'par', par, 'debug', debug);
            end
            
        case {'categorical', 'multinomial'}
            if isempty(c)
                ll = llCategoricalWarped(mu, f, type, 'loop', loop, 'par', par, 'debug', debug);
            else
                ll = llCategoricalPushed(mu, f, c, type, 'loop', loop, 'par', par, 'debug', debug);
            end
    end
    
end