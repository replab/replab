classdef IsomorphismComposition < replab.Isomorphism & replab.mrp.Composition

    methods

        function self = IsomorphismComposition(second, first, imageElementFun)
            self@replab.mrp.Composition(second, first, imageElementFun);
        end

    end

    methods % Implementations

        function s = preimageElement(t)
            s = first.preimageElement(second.preimageElement(t));
        end

        function t = imageElement(self, s)
            t = imageElement@replab.mrp.Composition(self, s);
        end
        
        % Obj

        function l = laws(self)
            l = laws@replab.Isomorphism(self);
        end
        
    end

    methods (Static)
        
        function m = lambda(source, target, preimageElementFun, imageElementFun, torusMap)
            m = lambda@replab.Isomorphism(source, target, preimageElementFun, imageElementFun, torusMap);
        end

    end
    
    methods (Access = protected) % Implementations

        function I = computeInverse(self)
            I = replab.mrp.compose(self.first.inverse, self.second.inverse);
        end

    end

end
