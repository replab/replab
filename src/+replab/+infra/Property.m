classdef Property < replab.Str

    properties
        name
        attributes
        docLines
    end

    methods

        function self = Property(name, attributes, docLines)
            self.name = name;
            self.attributes = attributes;
            self.docLines = docLines;
        end

    end

end
