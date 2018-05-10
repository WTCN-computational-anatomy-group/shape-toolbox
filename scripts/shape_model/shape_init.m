function [dat, model] = shape_init(dat, model, opt)
% Initialise all variables (that need it) 
% + lower bound stuff
    
    if opt.ui.verbose
        fprintf('%-20s%7s | %s |\n', 'VEM', 'Init', repmat('=',1,48));
    end

    model.lb = struct;

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
    
    if opt.ui.verbose
        fprintf('%-27s | ', 'Init Subjects #1');
    end
    
    % ---
    % Try to estimate max memory consumption
    optmem = opt.par.subjects.job.est_mem;
    maxmem = opt.par.subjects.job.mem;
    maxmem1 = (prod(opt.tpl.lat)*3*4*10)/(1024^2)+1.5; % 1.5G + 10 velocity fields
    maxmem1 = ceil(maxmem1 * 10)/10;
    opt.par.subjects.job.mem     = [num2str(maxmem1) 'G'];
    opt.par.subjects.job.est_mem = false;
    % ---
    
    [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
        'shape_init_subject', 'inplace', dat, model, opt);
    
    % ---
    opt.par.subjects.job.est_mem = optmem;
    opt.par.subjects.job.mem     = maxmem;
    % ---
    
    % ---------------------------------------------------------------------
    %    Template
    % ---------------------------------------------------------------------
    % These parameters depend on the above subject parameters
    % ---------------------------------------------------------------------
    
    if ~opt.tpl.provided && opt.f.N && opt.optimise.tpl.a
        if opt.ui.verbose
            fprintf('%-27s | \n', 'Init Template');
        end
        switch lower(opt.tpl.update)
            case 'map'
                model = aggregateTemplateGradHess(dat, model, opt);
                model = updateTemplate(model, opt);
                model = lbTemplate(model, opt);
            case 'ml'
                model = updateTemplateML(dat, model, opt);
        end
        model = updateTemplateDerivatives(model, opt);
    end
    
    % ---------------------------------------------------------------------
    %    Individual parameters (post-Template)
    % ---------------------------------------------------------------------
    
    if opt.ui.verbose
        fprintf('%-27s | ', 'Init Subjects #2');
    end
    
    % ---
    % Try to estimate max memory consumption
    optmem = opt.par.subjects.job.est_mem;
    maxmem = opt.par.subjects.job.mem;
    maxmem1 = (prod(opt.tpl.lat)*3*4*10)/(1024^2)+1.5; % 1.5G + 10 velocity fields
    maxmem1 = ceil(maxmem1 * 10)/10;
    opt.par.subjects.job.mem     = [num2str(maxmem1) 'G'];
    opt.par.subjects.job.est_mem = false;
    % ---
    
    [opt.par.subjects, dat] = distribute(opt.par.subjects, ...
        'shape_init_subject2', 'inplace', dat, model, opt);
    
    % ---
    opt.par.subjects.job.est_mem = optmem;
    opt.par.subjects.job.mem     = maxmem;
    % ---
    
    % ---------------------------------------------------------------------
    %    Aggregate
    % ---------------------------------------------------------------------
    
    if opt.optimise.q.q
        model = aggregateAffine(dat, model, opt);
    end
    if opt.optimise.z.z
        model = aggregateLatent(dat, model, opt);
    end
    if opt.optimise.v.v
        model = aggregateVelocity(dat, model, opt);
    end
    if opt.f.N
        model = aggregateMatching(dat, model, opt);
    end
    
    % ---------------------------------------------------------------------
    %    Update lower bound
    % ---------------------------------------------------------------------
    model = updateLowerBound(model);    % Accumulate lower bound parts
    shape_plot_all(model, opt);
    model.converged = false;
    
    % ---------------------------------------------------------------------
    %    Should decide what to do with that...
    % ---------------------------------------------------------------------
    model.q.active = true;
    model.v.active = true;
    model.pg.active = true;
end
