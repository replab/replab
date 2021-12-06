classdef SemidirectProductGroupType < replab.gen.StaticFiniteGroupType
% Describes the type of a semidirect product of finite groups
%
% For semidirect product groups, the `+replab.FiniteGroupType` machinery is overkill, but it still reduces code duplication.
%
% The nice isomorphism is computed using Proposition 2.1 from https://www.maths.usyd.edu.au/u/pubs/publist/preprints/2015/easdown-19.pdf
% (the published version does not contains the proof which interests us).

    properties (SetAccess = protected)
        H % (`+replab.FiniteGroup`): Group acting
        N % (`+replab.FiniteGroup`): Group acted upon
        phi % (`+replab.Action`): Action of H on N
        regularIsomorphism % (`+replab.FiniteIsomorphism`): Regular isomorphism of `.N`
        regularPhi % (`+replab.FiniteMorphism`): `.phi` as a morphism from `.H` to the symmetric group of degree ``N.order``
    end

    methods

        function self = SemidirectProductGroupType(phi)
            assert(isa(phi, 'replab.Action'));
            % note the code duplication with `.SemidirectProductGroup_compact` and `.SemidirectProductGroup_finite`
            H = phi.G;
            N = phi.P;
            assert(isa(N, 'replab.FiniteGroup'));
            assert(isa(H, 'replab.FiniteGroup'));
            self.identity = {H.identity N.identity};
            self.phi = phi;
            self.H = H;
            self.N = N;
            % compute the regular representation stuff
            permIso = N.orderPreservingPermutationIsomorphism;
            S = replab.perm.Set.fromPermutationGroup(permIso.target);
            o = double(N.order);
            % regular representation of N
            regularImages = cell(1, N.nGenerators);
            for i = 1:N.nGenerators
                img = zeros(1, o);
                gen = permIso.imageElement(N.generator(i));
                for j = 1:o
                    img(j) = S.find(gen(S.at(j)')');
                end
                regularImages{i} = img;
            end
            self.regularIsomorphism = N.isomorphismByImages(replab.S(o), 'preimages', N.generators, 'images', regularImages);
            % regular representation of phi
            phiImages = cell(1, H.nGenerators);
            for i = 1:H.nGenerators
                img = zeros(1, o);
                genH = H.generator(i);
                for j = 1:o
                    preimg_prm = S.at(j)';
                    preimg = permIso.preimageElement(preimg_prm);
                    img = phi.leftAction(genH, preimg);
                    img_prm = permIso.imageElement(img);
                    phiImage(j) = S.find(img_prm');
                end
                phiImages{i} = phiImage;
            end
            self.regularPhi = H.morphismByImages(replab.S(o), 'preimages', H.generators, 'images', phiImages);
            generators = horzcat(cellfun(@(h) {h, N.identity}, H.generators, 'uniform', 0), ...
                                 cellfun(@(n) {H.identity, n}, N.generators, 'uniform', 0));
            generatorNames = horzcat(H.generatorNames, N.generatorNames);
            if length(unique(generatorNames)) ~= length(generatorNames) % names are not unique
                generatorNames = [];
            end
            args = {'order', H.order*N.order, 'generatorNames', generatorNames};
            niceType = replab.PermutationGroupType.make(H.orderPreservingPermutationIsomorphism.target.domainSize + double(N.order));
            self.finishConstruction(generators, args, niceType);
        end

    end

% $$$     methods (Access = protected)
% $$$
% $$$         function iso = computeRegularPhi(self)
% $$$             o = double(self.N.order);
% $$$             nGenH = self.H.nGenerators;
% $$$             phiImages = cell(1, nGenH);
% $$$             s = self.N.elementsSequence;
% $$$             for i = 1:nGenH
% $$$                 phiImage = zeros(1, o);
% $$$                 genH = self.H.generator(i);
% $$$                 for j = 1:o
% $$$                     preimg = s.at(j);
% $$$                     img = self.phi.leftAction(genH, preimg);
% $$$                     phiImage(j) = s.find(img);
% $$$                 end
% $$$                 phiImages{i} = phiImage;
% $$$             end
% $$$             iso = self.H.morphismByImages(replab.S(o), 'preimages', self.H.generators, 'images', phiImages);
% $$$         end
% $$$
% $$$     end
% $$$
% $$$     methods
% $$$
% $$$         function iso = regularPhi(self)
% $$$         % Computes the action of `.phi` on the regular representation of `.N`
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.FiniteMorphism`: Morphism from `.H` to the symmetric group acting on ``N.order`` elements
% $$$             iso = self.cached('regularPhi', @() self.computeRegularPhi);
% $$$         end
% $$$
% $$$         function iso = regularIsomorphism(self)
% $$$         % Returns the regular isomorphism of `.N`
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.FiniteIsomorphism`: Isomorphism from `.N` to its regular permutation group
% $$$             iso = self.cached('regularIsomorphism', @() self.N.regularIsomorphism);
% $$$         end
% $$$
% $$$     end


    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end

        function g = sample(self)
            g = {self.H.sample, self.N.sample};
        end

        % Monoid

        function z = compose(self, x, y)
            xh = x{1};
            xn = x{2};
            yh = y{1};
            yn = y{2};
            % see `+replab.SemidirectProductGroup`
            yhinv = self.H.inverse(yh);
            zh = self.H.compose(xh, yh);
            zn = self.N.compose(self.phi.leftAction(yhinv, xn), yn);
            z = {zh zn};
        end

        % Group

        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = self.N.inverse(self.phi.leftAction(xh, xn));
            z = {zh zn};
        end

        % FiniteGroupType

        function l = isSameTypeAs(self, otherType)
            l = false; % guilty until proven innocent
            if ~isa(otherType, 'replab.prods.SemidirectProductGroupType')
                return
            end
            if ~(self.H == otherType.H) || ~(self.N == otherType.N)
                return
            end
            phi1 = self.phi;
            phi2 = otherType.phi;
            for i = 1:self.H.nGenerators
                h = self.H.generator(i);
                for j = 1:self.N.nGenerators
                    n = self.N.generator(j);
                    img1 = phi1.leftAction(h, n);
                    img2 = phi2.leftAction(h, n);
                    if ~self.N.eqv(img1, img2)
                        return
                    end
                end
            end
            l = true;
        end

        % StaticFiniteGroupType

        function t = imageElement(self, s)
            h = s{1};
            n = s{2};
            p1 = self.H.orderPreservingPermutationIsomorphism.imageElement(h);
            p2a = self.regularPhi.imageElement(h);
            p2b = self.regularIsomorphism.imageElement(n);
            p2 = p2a(p2b);
            t = [p1 p2+length(p1)];
        end

        function G = makeParentGroup(self, generators, nice, niceIsomorphism)
            G = replab.prods.SemidirectProductGroup_finite(self, generators, nice, niceIsomorphism);
        end

        function s = preimageElement(self, t)
            ds = self.H.orderPreservingPermutationIsomorphism.target.domainSize;
            p1 = t(1:ds);
            p2 = t((ds+1):end) - ds;
            h = self.H.orderPreservingPermutationIsomorphism.preimageElement(p1);
            p2a_inv = self.regularPhi.imageElement(self.H.inverse(h));
            p2b = p2a_inv(p2);
            n = self.regularIsomorphism.preimageElement(p2b);
            s = {h, n};
        end

    end

end
