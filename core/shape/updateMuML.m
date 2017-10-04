function mu = updateMuML(model, varargin)
% FORMAT mu = updateMuCategoricalML(model, pf1, ..., pfN, c1, ..., cN,
%                                   ('loop', loop), ('par', par))
% ** Required **
% model - Structure with field 'name' describing the generative model.
% pfi   - Image pushed in template space
% ci    - Pushed voxel count
% ** Keyword arguments **
% loop  - How to split processing: 'slice', 'none' or '' [auto]
% par   - Distribute compute [auto]
% ** Output **
% mu    - Updated template
%
% Closed form M-step update of the template.
% (Maximum likelihood update: no prior on Mu)

    switch lower(model.name)
        case {'normal', 'gaussian', 'l2'}
            mu = updateMuNormalML(varargin{:});
        case {'laplace', 'l1'}
            mu = updateMuLaplaceML(varargin{:});
        case {'bernoulli', 'binomial'}
            mu = updateMuBernoulliML(varargin{:});
        case {'categorical', 'multinomial'}
            mu = updateMuCategoricalML(varargin{:});
    end

end