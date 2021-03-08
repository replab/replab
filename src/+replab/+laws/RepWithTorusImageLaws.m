classdef RepWithTorusImageLaws < replab.Laws

    properties (SetAccess = protected)
        rep % (`+replab.Rep`): Representation
        mu % (`+replab.Morphism`): Morphism from maximal torus to ``rep.group``
        T % (`+replab.TorusGroup`): Maximal torus subgroup
    end

    methods

        function self = RepWithTorusImageLaws(rep)
            self.rep = rep;
            mu = rep.group.reconstruction;
            self.mu = mu;
            self.T = mu.source;
        end


        function law_torusImage_T(self, t)
            img1 = self.rep.image(self.mu.imageElement(t));
            [torusMap, torusInjection, torusProjection] = self.rep.torusImage;
            img2 = torusInjection * replab.TorusGroup.torusRepImage(torusMap * t) * torusProjection;
            if all(all(torusInjection == torusProjection'))
                c = 1;
            elseif issparse(torusInjection) && issparse(torusProjection)
                c = normest(torusProjection) * normest(torusProjection);
            else
                c = cond(full(torusInjection));
            end
            tol = 1e-14*c*sqrt(self.rep.dimension); % assumption: error on each phase is max 1e-15
            self.assertApproxEqual(img1, img2, tol);
        end

    end

end
