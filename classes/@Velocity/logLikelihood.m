function ll = logLikelihood(obj)
% FORMAT (ll) = obj.logLikelihood()
% ll - Model log-likelihood [if no argout: write to obj.ll]

    if nargout == 0 && obj.utd.ll
        ll = obj.ll;
        return
    end
    if obj.Debug, fprintf('* logLikelihood\n'); end;

    if nargout == 0
        obj.logLikelihoodMatching();
        obj.logLikelihoodLaplace();
        obj.ll = obj.llm + obj.lll;
        if obj.checkarray('z')
            obj.logLikelihoodPriorZ();
            obj.ll = obj.ll + obj.llz;
        end
        if obj.checkarray('r')
            obj.logLikelihoodPriorR();
            obj.ll = obj.ll + obj.llr;
        end
        obj.statusChanged('ll');
        obj.utd.ll = true;
    else
        ll = obj.logLikelihoodMatching() + obj.logLikelihoodLaplace();
        if obj.checkarray('z')
            ll = ll + obj.logLikelihoodPriorZ();
        end
        if obj.checkarray('r')
            ll = ll + obj.logLikelihoodPriorR();
        end
    end
end