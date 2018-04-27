function [dat, model] = shape_init(dat, model, opt)
% Initialise all variables (that need it) 
% + lower bound stuff
    

    % ---------------------------------------------------------------------
    %    Model parameters
    % ---------------------------------------------------------------------
    % These parameters can be provided or initialised from scratch.
    %  > For file arrays (w, a), this is specified by the [opt.provided]
    %    structure.
    %  > For precision parameters, they are always initialised from the
    %    provided prior value (opt.z.A0, opt.q.A0, opt.v.l0)
    % They can be considered parameters to optimise or be set fixed
    %  > This is specified by the [opt.optimise] structre
    % ---------------------------------------------------------------------
    
    model = shape_init_population(model, opt);
    
    % ---------------------------------------------------------------------
    %    Individual parameters (pre-Template)
    % ---------------------------------------------------------------------
    % These parameters cannot be provided and must be optimised.
    % ---------------------------------------------------------------------
    
    % Try to estimate max memory consumption
    optmem = opt.par.subjects.job.est_mem;
    maxmem = opt.par.subjects.job.mem;
    maxmem1 = (prod(opt.tpl.lat)*3*4*10)/(1024^2)+1.5; % 1.5G + 10 velocity fields
    maxmem1 = ceil(maxmem1 * 10)/10;
    opt.par.subjects.job.mem     = [num2str(maxmem1) 'G'];
    opt.par.subjects.job.est_mem = false;
    
    [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
        'shape_init_subject', 'iter', dat, model, opt);
    
    opt.par.subjects.job.est_mem = optmem;
    opt.par.subjects.job.mem     = maxmem;
    
    % ---------------------------------------------------------------------
    %    Template
    % ---------------------------------------------------------------------
    % These parameters depend on the above subject parameters
    % ---------------------------------------------------------------------
    
    if opt.f.N && opt.optimise.tpl.a
        switch lower(opt.tpl.update)
            case 'map'
                model = aggregateTemplateGradHess(dat, model, opt);
                model = updateTemplate(model, opt);
            case 'ml'
                model = updateTemplateML(dat, model, opt);
        end
    end
    
    % ---------------------------------------------------------------------
    %    Individual parameters (post-Template)
    % ---------------------------------------------------------------------
    
    % Try to estimate max memory consumption
    optmem = opt.par.subjects.job.est_mem;
    maxmem = opt.par.subjects.job.mem;
    maxmem1 = (prod(opt.tpl.lat)*3*4*10)/(1024^2)+1.5; % 1.5G + 10 velocity fields
    maxmem1 = ceil(maxmem1 * 10)/10;
    opt.par.subjects.job.mem     = [num2str(maxmem1) 'G'];
    opt.par.subjects.job.est_mem = false;
    
    [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
        'shape_init_subject2', 'iter', dat, model, opt);
    
    opt.par.subjects.job.est_mem = optmem;
    opt.par.subjects.job.mem     = maxmem;
    
    % ---------------------------------------------------------------------
    %    Aggregate
    % ---------------------------------------------------------------------
    
    model = aggregateAffine(dat, model, opt);
    model = aggregateLatent(dat, model, opt);
    model = aggregateVelocity(dat, model, opt);
    
    % ---------------------------------------------------------------------
    %    Update lower bound
    % ---------------------------------------------------------------------
    model = updateLowerBound(model);    % Accumulate lower bound parts
end
