classdef MethodGroup < replab.Str
% A group of methods in a class

    properties
        name % (charstring): Name of the group
        methodsInGroup % (cell(1,\*) of `.ClassElement`): Methods in the group
    end

    methods

        function self = MethodGroup(name, methodsInGroup)
            self.name = name;
            self.methodsInGroup = methodsInGroup;
        end

        function l = hasAccessibleMethods(self)
            l = any(cellfun(@(m) m.isAccessible, self.methodsInGroup));
        end

    end


end