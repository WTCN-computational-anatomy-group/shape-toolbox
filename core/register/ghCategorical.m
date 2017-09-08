function [g, h, htype] = ghCategorical(mu, f, c, varargin)
% FORMAT [g, (h, htype)] = ghCategorical(mu, f, c, (ga),
%                                   ('loop', loop), ('par', par))
% mu   - Reconstructed probability template
% f    - Observed image pushed in the template space.
% c    - Pushed voxel count.
% ga   - Spatial gradients of the log-probability template.
% loop - Specify how to split data processing
%        ('slice' or 'none' [default])
% par  - If true, parallelise processing [default: false]
% g    - Gradient
% h    - Hessian
% htype - Shape of the hessian ('diagonal', 'symtensor')
%
% Gradient & Hessian of the **negative** log-likelihood of the Categorical 
% matching term w.r.t. changes in the initial velocity.
%
% Let us remind that the matching term is of the form
% > E = sum c.F(f/c, mu)
% Hence,
% > dE = c.F'(f/c, mu) * dmu
% > d2E = tr(dmu) * c.F''(f/c, mu) * dmu
%
% If gmu is provided
% > Returns gradiant/hessian of the matching term w.r.t. changes in v.
% > size(g) = [nx ny nz 3]
% > size(h) = [nx ny nz 6]
% > htype   = 'symtensor'
% Else
% > Returns only gradients and hessian of the function F of f and mu
% > size(g) = [nx ny nz nc]
% > size(h) = [nx ny nz nc]
% > htype   = 'diagonal'

    % --- Parse inputs
    p = inputParser;
    p.FunctionName = 'ghCategorical';
    p.addRequired('mu',  @checkarray);
    p.addRequired('f',   @checkarray);
    p.addRequired('c',   @checkarray);
    p.addOptional('gmu', []);
    p.addParameter('loop',   '',    @ischar);
    p.addParameter('par',    false, @isscalar);
    p.addParameter('output', []);
    p.addParameter('debug',  false, @isscalar);
    p.parse(mu, f, c, varargin{:});
    gmu = p.Results.gmu;
    par = p.Results.par;
    loop = p.Results.loop;
    
    if p.Results.debug, fprintf('* ghCategorical\n'); end;
    
    % --- Optimise parallelisation and splitting schemes
    [par, loop] = autoParLoop(par, loop, isa(mu, 'file_array'), size(mu, 3));
    
    % --- Read dimensions
    dim  = [size(mu) 1 1];
    dlat = dim(1:3);
    nc   = dim(4);
    nvec = size(gmu, 5);
    
    % --- Prepare output
    output = p.Results.output;
    if ~iscell(output)
        output = {output};
    end
    if numel(output) < 2
        output = [output {[]}];
    end
    if isempty(gmu)
        g = prepareOnDisk(output{1}, [dlat nc]);
    else
        g = prepareOnDisk(output{1}, [dlat nvec]);
    end
    if nargout > 1
        if isempty(gmu)
            h = prepareOnDisk(output{2}, [dlat nc*(nc+1)/2]);
        else
            h = prepareOnDisk(output{2}, [dlat nvec*(nvec+1)/2]);
        end
    end
    
    % --- No loop
    if strcmpi(loop, 'none')
        if p.Results.debug, fprintf('   - No loop\n'); end;
        if isempty(gmu)
            if nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c);
            else
                g(:,:,:,:) = onMemory(mu, f, c);
            end
        else
            if nargout > 1
                [g(:,:,:,:), h(:,:,:,:)] = onMemory(mu, f, c, gmu);
            else
                g(:,:,:,:) = onMemory(mu, f, c, gmu);
            end
        end
    % --- Loop on slices
    elseif strcmpi(loop, 'slice')
        % - Parallelised
        if p.Results.debug
            if par > 0
                fprintf('   - Parallelise on slices\n'); 
            else
                fprintf('   - Serialise on slices\n'); 
            end
        end
        if isempty(gmu)
            if nargout > 1
                parfor (z=1:dlat(3), par)
                    [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z));
                    g(:,:,z,:) = g(:,:,z,:) + g1;
                    h(:,:,z,:) = h(:,:,z,:) + h1;
                end
             else
                parfor (z=1:dlat(3), par)
                    g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z));
                    g(:,:,z,:) = g(:,:,z,:) + g1;
                end
            end
        else
            if nargout > 1
                g = numeric(g);
                h = numeric(h);
                parfor (z=1:dlat(3), par)
                    [g1, h1] = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z), gmu(:,:,z,:,:));
                    g(:,:,z,:) = g(:,:,z,:) + g1;
                    h(:,:,z,:) = h(:,:,z,:) + h1;
                end
            else
                g = numeric(g);
                parfor (z=1:dlat(3), par)
                    g1 = onMemory(mu(:,:,z,:), f(:,:,z,:), c(:,:,z),  gmu(:,:,z,:,:));
                    g(:,:,z,:) = g(:,:,z,:) + g1;
                end
            end
        end
        
    end
    
    % --- Set hessian type
    htype = 'symtensor';
    
    % --- Write on disk
    if ~isempty(output{1})
        g = saveOnDisk(g, output{1}, 'name', 'g');
    end
    if nargout > 1 && ~isempty(output{2})
        h = saveOnDisk(h, output{2}, 'name', 'h');
    end
    
end

% --- Actual gradient and hessian computation
function [g, h] = onMemory(mu, f, c, gmu)
    
    if nargin < 4
        gmu = [];
    end
    lat = [size(mu) 1 1];
    lat = lat(1:3);
    nc  = size(mu, 4);
    
    mu = single(numeric(mu));
    f  = single(numeric(f));
    c  = single(numeric(c));
    
    g  = bsxfun(@times, c, mu) - f;
    if ~isempty(gmu)
        g = pointwise(gmu, g, 't');
    end
    if nargout > 1
        [ind, length] = symIndices(nc, 'n');
        h = zeros([lat length], 'like', g);
        % Diagonals of the tensors
        for k=1:nc
            h(:,:,:,ind(k,k)) = c .* (  mu(:,:,:,k) - mu(:,:,:,k).^2 );
        end
        % Upper/Lower part
        for l=1:nc
            for m=(l+1):nc
                h(:,:,:,ind(l,m)) = - c .* mu(:,:,:,l) .* mu(:,:,:,m);
            end
        end
        if ~isempty(gmu)
            nvec = size(gmu, 5);
            [indv, length] = symIndices(nvec, 'n');
            hh = h;
            h = zeros([lat length], 'like', h);
            for d=1:nvec
                for l=d:nvec
                    for k1=1:nc
                        for k2=1:nc
                            h(:,:,:,indv(d,l)) = h(:,:,:,indv(d,l)) ...
                                + hh(:,:,:,ind(k1,k2)) .* gmu(:,:,:,k1,d) .* gmu(:,:,:,k2,l);
                        end
                    end
                end
            end
        end
        
    end
    
    % Just in case
    g(~isfinite(g)) = 0;
    if nargout > 1
        h(~isfinite(h)) = 0;
    end
end