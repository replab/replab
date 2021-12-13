classdef FiniteComposition < replab.FiniteMorphism & replab.mrp.Composition

    methods

        function self = FiniteComposition(second, first, imageElementFun)
            self@replab.mrp.Composition(second, first, imageElementFun);
        end

        function s = preimageRepresentative(self, t)
            s = self.first.preimageRepresentative(self.second.preimageRepresentative(t));
        end

        function t = imageElement(self, s)
            if ~isempty(self.imageElementFun)
                f = self.imageElementFun;
                t = f(s);
                return
            end
            t = self.second.imageElement(self.first.imageElement(s));
        end

        % Str

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteMorphism(self);
        end
        
        % Obj

        function l = laws(self)
            l = laws@replab.FiniteMorphism(self);
        end
        
    end

    methods (Access = protected)

        function I = computeImage(self)
            I = self.second.imageGroup(self.first.image);
        end

        function K = computeKernel(self)
            K = self.first.preimageGroup(self.second.kernel);
        end

    end

end
