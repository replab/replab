classdef Isotypic < replab.Str
% Describes an isotypic component in the decomposition of a representation
    properties
        parent;
        copies;
        multiplicity;
        copyDimension;
        realType;
        group;
        field;
    end
    
    methods
        
        function self = Isotypic(parent, copies, realType)
            assert(length(copies) >= 1, 'Isotypic component cannot be empty');
            self.parent = parent;
            self.copies = copies;
            self.copyDimension = self.copies{1}.dimension;
            self.multiplicity = length(self.copies);
            self.realType = realType;
            self.group = parent.group;
            self.field = parent.field;
        end
        
        function s = headerStr(self)
            rt = self.realType;
            if isequal(rt, [])
                rt = '';
            end
            if self.multiplicity > 1
                s = sprintf('Isotypic component I(%d)x%s(%d)', self.multiplicity, rt, self.copyDimension);
            else
                s = sprintf('Isotypic component %s(%d)', rt, self.copyDimension);
            end
        end
    end
end
