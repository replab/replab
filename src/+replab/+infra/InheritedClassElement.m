classdef InheritedClassElement < replab.infra.ClassElement
% Describes an inherited class element, present by the virtue of a superclass
    
    methods

        function self = InheritedClassElement(codeBase, parentClass, name, kind, attributes)
            declaration = ''; % No declaration as it is inherited
            self = self@replab.infra.ClassElement(codeBase, parentClass, name, kind, declaration, attributes);
        end

    end

end
