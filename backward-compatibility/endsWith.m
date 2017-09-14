function sw = endsWith(s, pattern, varargin)
% [BACKWARD COMPATIBILITY]
%
% TF = endsWith(S, PATTERN)
% returns true if S ends with PATTERN. It is a lot less versatile than
% the MATLAB R2016 version: S and PATTERN must be character vectors (not
% cells).
%
% TF = endsWith(S, PATTERN, 'IgnoreCase', IGNORE)
% if IGNORE is true, ignores case when comparing S and PATTERN.

% - If bultin function exists, use it
if exist('startsWith', 'builtin')
    sw = startWith(s, pattern, varargin{:});
    return
end

% - Ignore case
if numel(varargin) >= 2 ...
        && strcmp(varargin{1}, 'IgnoreCase') ...
        && varargin{2} == true
    s = lower(s);
    pattern = lower(pattern);
end

% - Check endsWith
sw = length(s) >= length(pattern) && strcmp(s(length(s)-length(pattern)+1:length(s)), pattern);