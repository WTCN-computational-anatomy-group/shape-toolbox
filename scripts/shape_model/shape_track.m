function shape_track(shape_model)

    datenum = inf;
    while 1
        pause(10)
        stat = dir(shape_model);
        if stat.datenum ~= datenum
            [model, opt] = load(shape_model, 'model', 'opt');
            shape_plot_all(model, opt);
        end
    end


end