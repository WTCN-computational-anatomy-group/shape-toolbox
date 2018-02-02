function sw = startsWith(s, pattern, varargin)
% [BACKWARD COMPATIBILITY]
%
% TF = startsWith(S, PATTERN)
% returns true if S starts with PATTERN. It is a lot less versatile than
% the MATLAB R2016 version: S and PATTERN must be character vectors (not
% cells).
%
% TF = startsWith(S, PATTERN, 'IgnoreCase', IGNORE)
% if IGNORE is true, ignores case when comparing S and PATTERN.

% - If bultin function exists, use it
if exist('startsWith', 'builtin')
    sw = builtin('startsWith', s, pattern, varargin{:});
    return
end

% - Ignore case
if numel(varargin) >= 2 ...
        && strcmp(varargin{1}, 'IgnoreCase') ...
        && varargin{2} == true
    s = lower(s);
    pattern = lower(pattern);
end

% - Check startsWith
sw = length(s) >= length(pattern) && strcmp(s(1:length(pattern)), pattern);