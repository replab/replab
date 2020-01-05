classdef ConcreteClassElement < replab.infra.SourceElement & replab.infra.ClassElement
% An element of a class that has a declaration in its class

    methods

        function self = ConcreteClassElement(codeBase, package, parentClass, name, startLineNumber, kind, declaration, ...
                                             attributes, docLines, docLineNumbers)
            sourceIdentifier = parentClass.name;
            self = self@replab.infra.SourceElement(codeBase, package, sourceIdentifier, startLineNumber, ...
                                                   name, docLines, docLineNumbers);
            self = self@replab.infra.ClassElement(codeBase, parentClass, name, kind, declaration, attributes);
        end

        function [packagePath elementPath] = splitPath(self)
            packagePath = self.parentClass.package.packagePath;
            elementPath = {self.parentClass.name self.name};
        end

    end

end
