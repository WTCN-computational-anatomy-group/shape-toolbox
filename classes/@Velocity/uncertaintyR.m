function sr = uncertaintyR(obj)

    if nargin == 1 && obj.checkarray('sr')
        sr = obj.sr;
        return
    end
    
    if ~obj.checkarray('hr')
        if obj.Debug
            warning('Cannot compute R''s uncertainty. Missing arrays.')
        end
        return
    end
    if obj.Debug, fprintf('* uncertaintyR\n'); end;
    
    hr = obj.hr;
    
    sr = pointwise3(hr, 'i');
    obj.sr.dim = size(sr);
    obj.sr(:) = sr(:);
    obj.statusChanged('sr');
    sr = obj.sr;
    
end