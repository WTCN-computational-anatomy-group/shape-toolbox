function sr = uncertaintyR(obj)

    if ~obj.checkarray('hr')
        if obj.Debug
            warning('Cannot compute R''s uncertainty. Missing arrays.')
        end
        return
    end
    if obj.Debug, fprintf('* uncertaintyR\n'); end;
    
    hr = obj.hr;
    
    sr = pointwiseInv(hr, true);
    obj.sr.dim = size(sr);
    obj.sr(:) = sr(:);
    obj.statusChanged('sr');
    obj.utd.sr = true;
    sr = obj.sr;
    
end