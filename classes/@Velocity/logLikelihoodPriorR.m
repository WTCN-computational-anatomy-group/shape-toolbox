function ll = logLikelihoodPriorR(obj, r, fast)
% FORMAT ll = obj.logLikelihoodPriorZ((r), (fast))
% (r)    - Latent coordinates [default: obj.r]
% (fast) - Only compute r-dependant parts [default: false]
% ll     - Log-likelihood of the input coordinates
%
% Returns log p(R)
% 

    if nargin == 1 && obj.utd.llr
        ll = obj.llr;
        return
    end

    % -- Default arguments
    if nargin < 3
        fast = false;
        if nargin < 2
            r = obj.r;
        end
    end
    
    % --- Check all arrays are ready to be used
    if ~obj.checkarray(r)
        ll = 0;
        if obj.Debug
            warning('Cannot compute log-likelihood: missing arrays\n');
        end
        return
    end
    if obj.Debug, fprintf('* logLikelihoodPriorR\n'); end;
    
    % --- Load in memory
    r = single(numeric(r));
    dim = [size(r) 1 1];
    lat = dim(1:3);
    
    % --- Discard bad voxels (they should not exist but just in case)
    if ~fast, count = sumall(isfinite(r)); end;
    r(~isfinite(r)) = 0;
    
    % --- Compute log-likelihood
    m = spm_diffeo('vel2mom', r, [obj.VoxelSize obj.RegParam]);
    ll = r(:)' * m(:);
    
    % --- Add constants w.r.t. R
    if ~fast
        ll = ll + count * log(2 * pi);
        
        [~, ld] = spm_shoot_greens('kernel', lat, [obj.VoxelSize obj.RegParam]);
        ll = ll - 0.5*ld(1);
%         if obj.Debug
%             warning(['LL: discard 0.5 * det(L). ' ...
%                      'Don''t know how to compute it (with the kernel?).']);
%         end;
    end
    
    % --- Common factor
    ll = - 0.5 * ll;

    if nargout == 0
        obj.llr     = ll;
        obj.utd.llr = true;
        obj.statusChanged('llr');
    end
    
end
