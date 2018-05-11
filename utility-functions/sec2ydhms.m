function str_end_2 = sec2ydhms(time, compact)
% FORMAT  str = sec2ydhms(duration, compact)
% duration - duration in seconds
% compact  - Use compact format [false]
%
% Convert a duration in seconds (obtained from tic/toc) to a character
% representation in terms of years, dat, hour, minute, seconds.
%
% Default representation: 3 hours, 15 minutes, 1 second
% Compact representation: 3h 15m 1s

    if nargin < 2
        compact = false;
    end
    
    dur = duration(0,0,time);
    elapsed = floor(years(dur));
    dur = dur - years(elapsed(end));
    elapsed = [elapsed floor(days(dur))];
    dur = dur - days(elapsed(end));
    elapsed = [elapsed floor(hours(dur))];
    dur = dur - hours(elapsed(end));
    elapsed = [elapsed floor(minutes(dur))];
    dur = dur - minutes(elapsed(end));
    elapsed = [elapsed floor(seconds(dur))];
    if compact
        units = {'y' 'd' 'h' 'm' 's'};
        space = '';
    else
        units = {'year' 'day' 'hour' 'minute' 'second'};
        space = ' ';
    end
    str_end_2 = '';
    for i=1:numel(elapsed)
        if elapsed(i) > 0
            str_end_2 = [str_end_2 sprintf('%d%s%s', elapsed(i), space, units{i})];
            if ~compact && elapsed(i) > 1
                str_end_2 = [str_end_2 's'];
            end
            if ~compact && sum(elapsed(i+1:end)) > 0
                str_end_2 = [str_end_2 ', '];
            end
        end
    end
    if sum(elapsed) == 0
        str_end_2 = [str_end_2 sprintf('< 0%s%s', space, units{5})];
    end

end