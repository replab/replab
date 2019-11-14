classdef RowMajorShape < replab.tensor.Shape
    
    properties
        asColumnMajor % replab.tensor.Shape: Tensor shape as a column major object
    end
    
    methods
        
        function self = RowMajorShape(dimensions, group, isOrderColumnMajor)
            assert(~isOrderColumnMajor);
            n = length(dimensions);
            self = self@replab.tensor.Shape(dimensions, group, isOrderColumnMajor);
            acmDimensions = fliplr(dimensions);
            acmGroup = group.leftConjugateGroup(n:-1:1);
            self.asColumnMajor = replab.tensor.Shape.make(acmDimensions, acmGroup, true);
        end
        
        function sub = indToSub(self, ind)
            sub = self.asColumnMajor.indToSub(ind);
            sub = fliplr(sub);
        end
            
        function ind = subToInd(self, sub)
            ind = self.asColumnMajor.subToInd(fliplr(sub));
        end
        
        function orbit = subOrbit(self, sub)
            orbit = self.asColumnMajor.subOrbit(fliplr(sub));
            orbit = fliplr(orbit);
        end
        
    end
    
end
