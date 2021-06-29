classdef CommutingMorphismsMorphism < replab.Morphism
% A morphism whose image consists of the product of images of commuting morphisms

    properties (SetAccess = protected)
        morphisms % (cell(1,\*) of `+replab.Morphism`): Morphisms
    end

    methods

        function self = CommutingMorphismsMorphism(source, target, morphisms)
            self.source = source;
            self.target = target;
            self.morphisms = morphisms;
            if isempty(morphisms)
                sd = source.maximalTorusDimension;
                td = target.maximalTorusDimension;
                if ~isempty(sd) && ~isempty(td)
                    tm = zeros(td, sd);
                else
                    tm = [];
                end
            elseif all(cellfun(@(m) ~isempty(m.torusMap), morphisms))
                tm = morphisms{1}.torusMap;
                for i = 2:length(morphisms)
                    tm = tm + morphisms{i}.torusMap;
                end
            else
                tm = [];
            end
            self.torusMap = tm;
        end

        function t = imageElement(self, s)
            if isempty(self.morphisms)
                t = self.target.identity;
            else
                t = self.morphisms{1}.imageElement(s);
                for i = 2:length(self.morphisms)
                    t = self.target.compose(t, self.morphisms{i}.imageElement(s));
                end
            end
        end

    end

end
