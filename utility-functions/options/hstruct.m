classdef hstruct < handle
	properties
		data
	end
	
	methods
		function obj = hstruct(data)
            if nargin < 1
                data = struct;
            end
			obj.data = data;
		end
	end
end