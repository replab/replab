classdef ClassElement < replab.Str
    
    properties
        name % charstring: Method or property identifier
        attributes % struct: Attributes from the ``methods``/``properties`` block
        doc % `.Doc`: Documentation
        packageNameParts % row cell vector of charstring: parent package identifier
        className % charstring: parent class name
        kind % {'property' 'method'}: Human readable description of type
        lineNumber % integer: First line of the element declaration 
    end
    
    methods
        
        function str = fieldNameEncoding(self)
            str = strjoin(horzcat(self.packageNameParts, {self.className self.name}), '_');
        end
        
        function str = fullName(self)
            str = strjoin(horzcat(self.packageNameParts, {self.className self.name}), '.');
        end

        function str = classFullName(self)
            str = strjoin(horzcat(self.packageNameParts, {self.className}), '.');
        end
        
        function str = sphinxFullName(self)
            pp = strjoin(cellfun(@(x) ['+' x], self.packageNameParts, 'uniform', 0), '.');
            str = [pp '.' self.className '.' self.name];
        end
        
    end
            
end
