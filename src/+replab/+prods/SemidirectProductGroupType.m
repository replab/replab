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
    end

    methods

        function self = SemidirectProductGroupType(phi)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            assert(isa(N, 'replab.FiniteGroup'));
            assert(isa(H, 'replab.FiniteGroup'));
            % note the code duplication with `.SemidirectProductGroup_compact` and `.SemidirectProductGroup_finite`
            self.identity = {H.identity N.identity};
            self.phi = phi;
            self.H = H;
            self.N = N;
        end

    end

    methods (Access = protected)

        function iso = computePhiMorphism(self)
            o = double(self.N.order);
            nGenH = self.H.nGenerators;
            phiImages = cell(1, nGenH);
            s = self.N.elementSequence;
            for i = 1:nGenH
                phiImage = zeros(1, o);
                genH = self.H.generator(i);
                for j = 1:o
                    preimg = s.at(j);
                    img = self.phi.leftAction(genH, preimg);
                    phiImage(j) = s.find(img);
                end
                phiImages{i} = phiImage;
            end
            iso = self.H.morphismByImages(replab.S(o), 'preimages', self.H.generators, 'images', phiImages);
        end

    end

    methods

        function iso = phiMorphism(self)
            iso = self.cached('phiMorphism', @() self.computePhiMorphism);
        end

    end


    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end

        function g = sample(self)
            g = {self.H.sample, self.N.sample};
        end

        % TotalOrder

        function c = compare(self, x, y)
            error('Abstract');
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
                    n = self.n.generator(j);
                    img1 = phi1.leftAction(h, n);
                    img2 = phi2.leftAction(h, n);
                    if ~self.N.eqv(img1, img2)
                        return
                    end
                end
            end
            l = true;
        end

    end

end
