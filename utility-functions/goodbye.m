function goodbye(global_start, name)
    
    global_end = toc(global_start);
    fprintf('\n');
    fprintf([' ' repmat('-',1,78) ' \n']);
    str_end_1 = sprintf('%20s || %10s ended.', name, datestr(now));
    fprintf(['| ' str_end_1 repmat(' ', 1, 80-3-length(str_end_1)) '|\n']);
    str_end_2 = sprintf('%20s || ', 'Elapsed time');
    % Convert to units
    dur = duration(0,0,global_end);
    elapsed = floor(years(dur));
    dur = dur - years(elapsed(end));
    elapsed = [elapsed floor(days(dur))];
    dur = dur - days(elapsed(end));
    elapsed = [elapsed floor(hours(dur))];
    dur = dur - hours(elapsed(end));
    elapsed = [elapsed floor(minutes(dur))];
    dur = dur - minutes(elapsed(end));
    elapsed = [elapsed floor(seconds(dur))];
    units   = {'year' 'day' 'hour' 'minute' 'second'};
    for i=1:numel(elapsed)
        if elapsed(i) > 0
            str_end_2 = [str_end_2 sprintf('%d %s', elapsed(i), units{i})];
            if elapsed(i) > 1
                str_end_2 = [str_end_2 's'];
            end
            if sum(elapsed(i+1:end)) > 0
                str_end_2 = [str_end_2 ', '];
            end
        end
    end
    if sum(elapsed) == 0
        str_end_2 = [str_end_2 '< 0 second'];
    end
    fprintf(['| ' str_end_2 repmat(' ', 1, 80-3-length(str_end_2)) '|\n']);
    str_goodbye = sprintf('%77s', 'Goodbye! ');
    fprintf(['| ' str_goodbye repmat(' ', 1, 80-3-length(str_goodbye)) '|\n']);
    fprintf([' ' repmat('-',1,78) ' \n\n']);
    diary off
    
end