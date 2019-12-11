classdef Property < replab.Str

    properties
        name % charstring: Property identifier
        attributes % struct: Attributes from the ``properties`` block
        docLines % row cell vector of charstring: Documentation comment lines
    end

    methods

        function self = Property(name, attributes, docLines)
            self.name = name;
            self.attributes = attributes;
            self.docLines = docLines;
        end
        
        function str = headerStr(self)
            keywords = {};
            if isfield(self.attributes, 'Access') && ismember(self.attributes.Access, {'private', 'protected'})
                keywords{1,end+1} = self.attributes.Access;
            end
            keywords{1,end+1} = 'property';
            str = [self.name ' (' strjoin(keywords, ' ') ')'];
        end

    end

end
