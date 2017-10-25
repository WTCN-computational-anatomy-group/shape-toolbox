function mu = updateMuML(model, dat, varargin)
% FORMAT mu = updateMuCategoricalML(model, dat, ('fwhm', fwhm),
%                                   ('loop', loop), ('par', par))
% ** Required **
% model - Structure with field 'name' describing the generative model.
% pfi   - Image pushed in template space
% ci    - Pushed voxel count
% ** Keyword arguments **
% fwhm  - Smoothing kernel used as pseudo prior [do not use]
% loop  - How to split processing: 'slice', 'none' or '' [auto]
% par   - Distribute compute [auto]
% ** Output **
% mu    - Updated template
%
% Closed form M-step update of the template.
% (Maximum likelihood update: no prior on Mu)

    switch lower(model.name)
        case {'normal', 'gaussian', 'l2'}
            if isfield(dat, 'sigma2')
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ...
                                      's', dat.sigma2, varargin{:});
            else
                mu = updateMuNormalML('f', dat.pf, 'c', dat.c, ...
                                      varargin{:});
            end
        case {'laplace', 'l1'}
            if isfield(dat, 'b')
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ...
                                       'b', dat.b, varargin{:});
            else
                mu = updateMuLaplaceML('f', dat.pf, 'c', dat.c, ... 
                                       varargin{:});
            end
        case {'bernoulli', 'binomial'}
            mu = updateMuBernoulliML(dat.pf, dat.c, varargin{:});
        case {'categorical', 'multinomial'}
            mu = updateMuCategoricalML(dat.pf, dat.c, varargin{:});
    end

end