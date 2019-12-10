classdef Method < replab.Str
    
    properties
        name
        attributes
        declaration
        docLines
    end
    
    methods
        
        function self = Method(name, attributes, declaration, docLines)
            self.name = name;
            self.attributes = attributes;
            self.declaration = declaration;
            self.docLines = docLines;
        end
        
    end
    
end
