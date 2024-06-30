function argstruct = ERATO_testcase(varargs)
%
% 
%   USAGE: argstruct = setargs(varargs)
% __________________________________________________________________________
% OUTPUT
% 
% 	ARGSTRUCT
%    structure containing the final argument values
% __________________________________________________________________________
% INPUTS
% 	VARARGS [optional]     
%     cell array of user-specified "'Name', value" pairs for one or more of
%     the variables with default values. 
% __________________________________________________________________________
% USAGE EXAMPLE 
%  
%     argstruct   = ERATO_matlab(varargin)

%% interpretation of the varargs
if nargin < 1, varargs = []; end
defaultargs = {'E',2;'a',1/3;'q0',0.3;'n',2;'n_s',14};
if ~isempty(varargs)
    if mod(length(varargs), 2)
        error('Optional inputs must be entered as "''Name'', Value" pairs, e.g., myfunction(''arg1'', val1, ''arg2'', val2)'); 
    end
    arg = reshape(varargs, 2, length(varargs)/2)';
    for i = 1:size(arg,1)
       idx = strncmpi(defaultargs(:,1), arg{i,1}, length(arg{i,1}));
       if sum(idx) > 1
           error(['Input "%s" matches multiple valid inputs:' repmat('  %s', 1, sum(idx))], arg{i,1}, defaultargs{idx, 1});
       elseif ~any(idx)
           error('Input "%s" does not match a valid input.', arg{i,1});
       else
           defaultargs{idx,2} = arg{i,2};
       end
    end
end
for i = 1:size(defaultargs,1)
    assignin('caller', defaultargs{i,1}, defaultargs{i,2})
end

end


