function methodList = methodList(mc)
% Returns the list of methods for the given metaclass
%
% Args:
%   mc (metaclass): Metaclass whose methods to enumerate
%
% Returns:
%   cell(1,\*) of method descriptions: The metaclass methods
    methodList = mc.MethodList;
    methodList = methodList(:).';
    if ~isa(methodList, 'cell')
        methodList = num2cell(methodList);
    end
end
