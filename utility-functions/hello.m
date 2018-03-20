function cleanupObj = hello(name)
    global_start = tic;
    fprintf([' ' repmat('-',1,78) ' \n']);
    str_hello = 'Hello!';
    fprintf(['| ' str_hello repmat(' ', 1, 80-3-length(str_hello)) '|\n']);
    str_started = sprintf('%20s || %10s started...', name, datestr(now));
    fprintf(['| ' str_started repmat(' ', 1, 80-3-length(str_started)) '|\n']);
    fprintf([' ' repmat('-',1,78) ' \n\n']);
    cleanupObj = onCleanup(@() goodbye(global_start, name));
end