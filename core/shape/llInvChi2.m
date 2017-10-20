function ll = llInvChi2(x, k, fast)
    
    if nargin < 3
        fast = false;
    else
        fast =strcmpi(fast, 'fast');
    end

    ll = - (k/2+1) .* log(x) - 1./(2*x);
    if ~fast
        ll = ll - k./2.*log(2) - gammaln(k/2);
    end

end