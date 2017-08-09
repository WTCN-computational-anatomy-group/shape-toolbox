function z = initZ(obj, method, varargin)
% FORMAT (z) = obj.initZ( (method, (args)) )
% (method) - 'zeros' [default], 'ones', 'normal', 'rand', 'discard'
% (args)   - if method is 'normal': (args[1]) - variance [default: 1]
%                                   (args[2]) - mean [default: 0]
% (z)      - Output vector [if no argout: write to obj.z]
%
% Initialise latent coordinates

    if obj.Debug, fprintf('* initZ\n'); end;
    
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
        obj.z.dim   = [0 0];
%         obj.hz.dim  = [0 0];
%         obj.sz.dim  = [0 0];
%         obj.swz.dim = [0 0];
        obj.statusChanged('hz', 'sz', 'swz');
    else
        % --- Latent dimension
        k = size(obj.w, 5);

        % --- Init
        switch method
            case 'zeros'
                z = zeros(k, 1, 'single');
            case 'ones'
                z = ones(k, 1, 'single');
            case 'normal'
                z = single(normrnd(mean, sqrt(var), k, 1));
            case 'rand'
                z = single(randn(k, 1));
        end

        % --- Write on disk
        if nargout == 0
            obj.z.dim = size(z);
            obj.z(:)  = z(:);
        end
    end
    if nargout == 0
        obj.statusChanged('z');
        obj.utd.z = true;
    end
        
    
end