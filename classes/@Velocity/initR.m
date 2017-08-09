function r = initR(obj, method, varargin)
% FORMAT (r) = obj.initR( (method, (args)) )
% (method) - 'zeros' [default], 'ones', 'normal', 'rand'
% (args)   - if method is 'normal': (args[1]) - variance [default: 1]
%                                   (args[2]) - mean [default: 0]
% (r)      - Output vector [if no argout: write to obj.r]
% obj      - Velocity object. This function uses property 'mu' (to read
%            dimensions info)
%
% Initialise residual field.


    obj.reconstructTemplate();
    if obj.Debug, fprintf('* initR\n'); end;
    
    % --- Default arguments
    if nargin < 2
        method = 'zeros';
    end
    switch method
        case 'normal'
            var  = 1;
            mean = 0;
            if numel(varargin) >= 1
                var = varargin{1};
                if numel(varargin) >= 2
                    mean = varargin{2};
                end
            end
    end
    
    if strcmpi(method, 'discard')
        % --- Discard case
        obj.r.dim  = [0 0];
%         obj.hr.dim = [0 0];
%         obj.sr.dim = [0 0];
        obj.statusChanged('hr', 'sr');
    else
        % --- Latent dimension
        dim = [size(obj.mu) 1];
        dim = [dim(1:3) 3];

        % --- Init
        switch lower(method)
            case 'zeros'
                r = zeros(dim, 'single');
            case 'ones'
                r = ones(dim, 'single');
            case 'normal'
                r = single(normrnd(mean, sqrt(var), dim));
            case 'rand'
                r = single(randn(dim));
        end
        % --- Write on disk
        if nargout == 0
            obj.r.dim = size(r);
            obj.r(:)  = r(:);
        end
    end
    if nargout == 0
        obj.statusChanged('r');
        obj.utd.r = true;
    end
    
end