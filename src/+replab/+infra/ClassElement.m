classdef ClassElement < replab.infra.Element
% Describes a class element, which can be either a property or a method, concrete or inherited
%
% Although it is not immediately evident, the class `.ClassElement` is abstract.

    properties
        parentClass % `.Class`: Parent class
        kind % {'method', 'property'}: Kind of the class member
        declaration % charstring: Declaration line for method, if available / ``[]`` for property
        attributes % struct: Attributes from the ``methods``/``properties`` block
    end

    methods

        function self = ClassElement(codeBase, parentClass, name, kind, declaration, attributes)
            self = self@replab.infra.Element(codeBase, name);
            self.parentClass = parentClass;
            self.kind = kind;
            self.declaration = declaration;
            self.attributes = attributes;
        end

        function [packagePath elementPath] = splitPath(self)
            packagePath = self.parentClass.package.packagePath;
            elementPath = {self.parentClass.name self.name};
        end

        function c = childrenNames(self)
            c = {};
        end

        function d = declarations(self)
            d = replab.infra.Declarations(self.codeBase, self);
        end

        function b = isAccessible(self)
            b = ~isfield(self.attributes, 'Access') || isequal(self.attributes.Access, 'public');
        end

        function b = isMethod(self)
            b = isequal(self.kind, 'method');
        end

        function b = isProperty(self)
            b = isequal(self.kind, 'property');
        end

    end

end
