function q = initAffine(obj, method)
% FORMAT (z) = obj.initZ( (method, (args)) )
% (method) - 'zeros' [default], 'discard'
% (q)      - Output vector [if no argout: write to obj.q]
%
% Initialise latent coordinates

    if obj.Debug, fprintf('* initQ\n'); end;
    
    % --- Default arguments
    if nargin < 2
        method = 'zeros';
    end
    
    if strcmpi(method, 'discard')
        q = [];
    else
        Q = size(obj.AffineBasis, 3);
        q = zeros(Q, 1);
    end
    
    if nargout == 0
        obj.q = q;
        obj.statusChanged('q');
    end
        
    
end