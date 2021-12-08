classdef FiniteIsomorphismComposition < replab.FiniteIsomorphism & replab.mrp.FiniteComposition & replab.mrp.IsomorphismComposition

    methods

        function self = FiniteIsomorphismComposition(second, first, imageElementFun)
            self@replab.mrp.FiniteComposition(second, first, imageElementFun);
            self@replab.mrp.IsomorphismComposition(second, first, imageElementFun);
        end

        function s = preimageElement(self, t)
            s = first.preimageElement(second.preimageElement(t));
        end

        function s = preimageRepresentative(self, t)
            s = preimageRepresentative@replab.mrp.FiniteComposition(self, t);
        end

        function l = preservesTypeOrder(self)
            l = self.first.preservesTypeOrder && self.second.preservesTypeOrder;
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            t = imageElement@replab.mrp.FiniteComposition(self, s);
        end

        % Str

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteIsomorphism(self);
        end

        % FiniteMorphism

        function T = imageGroup(self, S)
            T = imageGroup@replab.FiniteIsomorphism(self, S);
        end

        function S = preimageGroup(self, T)
            S = preimageGroup@replab.FiniteIsomorphism(self, T);
        end

        function m = restrictedSource(self, newSource)
            m = restrictedSource@replab.FiniteIsomorphism(self, newSource);
        end

        % Obj

        function l = laws(self)
            l = laws@replab.FiniteIsomorphism(self);
        end
        
    end
    
    methods (Static)
        
        function m = lambda(source, target, preimageElementFun, imageElementFun, torusMap)
            m = lambda@replab.mrp.IsomorphismComposition(source, target, preimageElementFun, imageElementFun, torusMap);
        end

        function m = identity(group)
            m = identity@replab.FiniteIsomorphism(group);
        end

    end
    
    methods (Access = protected)

        function I = computeInverse(self)
            I = replab.mrp.compose(self.first.inverse, self.second.inverse);
        end

        function K = computeKernel(self)
            K = computeKernel@replab.mrp.FiniteComposition(self);
        end

        function I = computeImage(self)
            I = self.target;
        end

    end

end
