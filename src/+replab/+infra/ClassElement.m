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

        function a = argumentString(self)
        % Returns the part of the method declaration that describes the method arguments
        %
        % Returns the empty string if this class element is a property, or if the method
        % is declared without arguments
        %
        % Returns:
        %   charstring: Argument declaration including parentheses
            a = '';
            if ~isempty(self.declaration)
                token = regexp(self.declaration, '(\([^\(]*\))', 'tokens', 'once');
                if isempty(token)
                    return
                end
                if iscell(token)
                    token = token{1};
                end
                a = token;
            end
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

        function decl = contextualizeDeclaration(self, fullIdentifier)
        % Returns the declaration, expressed with the given full name
            decl = strrep(self.declaration, self.name, fullIdentifier);
            switch self.kind
                case 'method'
                    if isequal(lower(decl(1:9)), 'function ')
                        decl = ['Method ', decl(10:end)];
                    end
            end
        end

    end

end
