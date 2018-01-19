function [dat, model] = translatePath(dat, model, onclient, onserver, serversep, clientsep)
% FORMAT [dat, model] = translatePath(dat, model, 
%                                     onclient,  onserver, 
%                                     serversep, clientsep)
% WARNING
% -------
% File separator conversion is extremely basic right now.
% It performs a blind substitution '/' <-> '\'
% When paths contain escaped characters (such as '\\'), this causes
% havoc.

    if nargin < 6
        clientsep = filesep;
        if nargin < 5
            serversep = filesep;
        end
    end

    datfields = fieldnames(dat);
    for n=1:numel(dat)
        for i=1:numel(datfields)
            f = datfields{i};
            if isa(dat(n).(f), 'file_array')
                dat(n).(f).fname = strrep(dat(n).(f).fname, onclient, onserver);
                if clientsep ~= serversep
                    dat(n).(f).fname = strrep(dat(n).(f).fname, clientsep, serversep);
                end
            end
        end
    end
    
   if ~isempty(model)
       modelfields = fieldnames(model);
       for i=1:numel(modelfields)
           f = modelfields{i};
           if isa(model.(f), 'file_array')
               model.(f).fname = strrep(model.(f).fname, onclient, onserver);
                if clientsep ~= serversep
                    model.(f).fname = strrep(model.(f).fname, clientsep, serversep);
                end
           end
       end
   end
   
end