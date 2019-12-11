classdef Method < replab.Str
    
    properties
        name
        attributes
        declaration
        docLines
        isAbstract
    end
    
    methods
        
        function self = Method(name, attributes, declaration, docLines, isAbstract)
            self.name = name;
            self.attributes = attributes;
            self.declaration = declaration;
            self.docLines = docLines;
            self.isAbstract = isAbstract;
        end
        
        function str = headerStr(self)
            keywords = {};
            if self.isAbstract
                keywords{1,end+1} = 'abstract';
            end
            if isfield(self.attributes, 'Access') && ismember(self.attributes.Access, {'private', 'protected'})
                keywords{1,end+1} = self.attributes.Access;
            end
            if isfield(self.attributes, 'Static') && self.attributes.Static
                keywords{1,end+1} = 'static';
            end
            keywords{1,end+1} = 'method';
            str = [self.name ' (' strjoin(keywords, ' ') ')'];
        end

    end
    
end
