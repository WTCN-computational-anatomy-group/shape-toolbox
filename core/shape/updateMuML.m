function mu = updateMuML(model, dat, varargin)
%__________________________________________________________________________
%
% Closed form M-step update of the template.
% (Maximum likelihood update: no prior on Mu)
%
% -------------------------------------------------------------------------
%
% FORMAT mu = updateMuML(model, dat, ...)
%
% REQUIRED
% --------
% model - Data model. Field 'name' must be one of
%         'normal', 'laplace', 'bionomial', 'categorical'
% dat   - Data structure
%
% KEYWORD ARGUMENTS
% -----------------
% lat  - Template lattice [temporarily REQUIRED]
% fwhm  - Smoothing kernel used as pseudo prior [do not use]
% loop  - How to split processing: 'slice', 'none' or '' [auto]
% par   - Distribute compute [auto]
%
% OUTPUT
% ------
% mu    - Updated template
%
%__________________________________________________________________________

    % Deal with both old and new factoring
    if ~isfield(dat, 'pf') && isstruct(dat(1).f)
        tmp = dat;
        dat = struct('pf', cell(numel(dat),1));
        for n=1:numel(dat)
            dat(n).pf = tmp(n).f.pf;
        end
        [dat.c] = deal([]);
        for n=1:numel(dat)
            dat(n).c = tmp(n).f.c;
        end
        if isfield(tmp(1).f, 'bb')
            [dat.bb] = deal([]);
            for n=1:numel(dat)
                dat(n).bb = tmp(n).f.bb;
            end
        end
        if isfield(tmp(1).f, 'sigma2')
            [dat.sigma2] = deal([]);
            for n=1:numel(dat)
                dat(n).sigma2 = tmp(n).f.sigma2;
            end
        end
        clear tpm
    end

    switch lower(model.name)
        case {'normal', 'gaussian', 'l2'}
            if isfield(dat, 'sigma2') && isfield(dat, 'bb')
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ...
                                      's', dat.sigma2, 'bb', dat.bb, ...  
                                      varargin{:});
            elseif isfield(dat, 'sigma2')
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ...
                                      's', dat.sigma2, ... 
                                      varargin{:});
            elseif isfield(dat, 'bb')
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ...
                                      'bb', dat.bb, ... 
                                      varargin{:});
            else
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ... 
                                      varargin{:});
            end
        case {'laplace', 'l1'}
            if isfield(dat, 'b') && isfield(dat, 'bb')
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ...
                                       'b', dat.b, 'bb', dat.bb, ... 
                                       varargin{:});
            elseif isfield(dat, 'b')
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ...
                                       'b', dat.b, ... 
                                       varargin{:});
            elseif isfield(dat, 'bb')
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ...
                                       'bb', dat.bb, ... 
                                       varargin{:});
            else
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ...  
                                       varargin{:});
            end
        case {'bernoulli', 'binomial'}
            if isfield(dat, 'bb')
                mu = updateMuBernoulliML('f', dat.pf, 'c', dat.c, ...
                                         'bb', dat.bb, ...
                                         varargin{:});
            else
                mu = updateMuBernoulliML('f', dat.pf, 'c', dat.c, ... 
                                         varargin{:});
            end
        case {'categorical', 'multinomial'}
            if isfield(dat, 'bb')
                mu = updateMuCategoricalML('f', dat.pf, 'c', dat.c, ...
                                           'bb', dat.bb, ...
                                            varargin{:});
            else
                mu = updateMuCategoricalML('f', dat.pf, 'c', dat.c, ... 
                                           varargin{:});
            end
    end

end