classdef ClassElement < replab.Str
    
    properties
        name % charstring: Method or property identifier
        attributes % struct: Attributes from the ``methods``/``properties`` block
        doc % `.Doc`: Documentation
        packageNameParts % row cell vector of charstring: parent package identifier
        className % charstring: parent class name
    end
    
    methods
        
        function str = fieldNameEncoding(self)
            str = strjoin(horzcat(self.packageNameParts, {self.className self.name}), '_');
        end
        
    end
            
end
