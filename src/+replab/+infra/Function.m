classdef Function < replab.infra.SourceElement
% Describes a MATLAB function
    
    properties
        declaration % charstring: Function declaration line
    end
    
    methods
        
        function self = Function(codeBase, package, functionData)
            self = self@replab.infra.SourceElement(codeBase, package, functionData.name, functionData.declarationLineNumber, ...
                                                   functionData.name, functionData.docLines, functionData.docLineNumbers);
            self.declaration = functionData.declaration;
        end
        
        function str = name(self)
            str = self.sourceIdentifier;
        end
        
        function p = elementPath(self)
            p = {self.name};
        end
               
        function c = childrenNames(self)
            c = {};
        end

        function decl = fullDeclaration(self)
        % Returns the declaration, expressed with the full identifier
            decl = strrep(self.declaration, self.name, self.fullIdentifier);
        end
        
    end
    
end
