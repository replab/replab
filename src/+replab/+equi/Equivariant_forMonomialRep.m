classdef Equivariant_forMonomialRep < replab.Equivariant
% Equivariant space of monomial representations of a finite group
%
% Example:
%   >>> g1 = [1 -2 -3 -4 -5];
%   >>> g2 = [1 4 5 2 3];
%   >>> g3 = [1 3 2 4 -5];
%   >>> G = replab.SignedPermutationGroup.of(g1,g2,g3);
%   >>> rep = G.naturalRep;
%   >>> isa(rep.commutant, 'replab.equi.Equivariant_forMonomialRep')
%       1

    properties (SetAccess = protected)
        phaseOrder % (integer): Phase order
        imagesR % (cell(1,\*) of generalized permutations): Generalized permutations acting on the row space (non-empty)
        imagesC % (cell(1,\*) of generalized permutations): Generalized permutatioms acting on the column space (non-empty)
    end

    methods

        function self = Equivariant_forMonomialRep(repR, repC, special, phaseOrder, imagesR, imagesC)
            self@replab.Equivariant(repR, repC, special);
            self.phaseOrder = phaseOrder;
            self.imagesR = imagesR;
            self.imagesC = imagesC;
        end

        function P = phasedMatrixPartition(self)
        % Returns the phased matrix partition corresponding to this equivariant space
        %
        % Returns:
        %   `.PhasedMatrixPartition`: The corresponding phased matrix partition
            P = self.cached('phasedMatrixPartition', @() self.computePhasedMatrixPartition);
        end

    end

    methods (Access = protected)

        function P = computePhasedMatrixPartition(self)
            P = replab.equi.PhasedMatrixPartition.fromGeneralizedPermutations(self.phaseOrder, self.imagesR, self.imagesC);
        end

    end

    methods % Implementations

        % Equivariant

        function b = isExact(self)
            b = true;
        end

    end

    methods (Static)

        function E = make(repR, repC, special)
            E = [];
            group = repR.group;
            preimages = group.generators;
            if isempty(preimages)
                preimages = {group.identity};
            end
            [GR, imagesR] = replab.rep.monomialImages(repR, preimages);
            if isempty(GR)
                return
            end
            [GC, imagesC] = replab.rep.monomialImages(repC, preimages);
            if isempty(GC)
                return
            end
            po = lcm(GR.m, GC.m);
            fR = po/GR.m;
            fC = po/GC.m;
            imagesR = cellfun(@(im) [im(1,:); fR*im(2,:)], imagesR, 'uniform', 0);
            imagesC = cellfun(@(im) [im(1,:); fC*im(2,:)], imagesC, 'uniform', 0);
            E = replab.equi.Equivariant_forMonomialRep(repR, repC, special, po, imagesR, imagesC);
        end

    end


    methods (Access = protected) % Implementations

        % Equivariant

        function X1 = project_exact(self, X)
            P = self.phasedMatrixPartition;
            nR = self.nR;
            nC = self.nC;
            X1 = replab.cyclotomic.zeros(nR, nC);
            po = P.phaseOrder;
            phases = replab.cyclotomic.zeros(1, po);
            for i = 1:po
                phases(i) = replab.cyclotomic.E(po)^(i-1);
            end
            for i = 1:P.nBlocks
                blk = P.blocks{i};
                n = size(blk, 2);
                ind = blk(1,:) + nR*(blk(2,:)-1);
                ph = P.phase(ind);
                coeffs = X(ind).*conj(phases(ph+1)); % multiply by conjugate phases
                s = sum(coeffs)/n;
                X1(ind) = s * phases(ph+1); % set the coefficients
            end
        end

        function [X1, eX1] = project_double_sparse(self, X)
            P = self.phasedMatrixPartition;
            nR = self.nR;
            nC = self.nC;
            X1 = zeros(nR, nC);
            E = zeros(nR, nC);
            po = P.phaseOrder;
            % compute roots of unity
            phases = [1 exp(2i*(1:po-1)*pi/po)];
            % one round of Newton iteration for additional precision
            phases = phases - phases .* (1 - phases.^(-po))/po;
            % force -1 when it's there
            if mod(po, 2) == 0
                phases(po/2+1) = -1;
            end
            for i = 1:P.nBlocks
                blk = P.blocks{i};
                n = size(blk, 2);
                ind = blk(1,:) + nR*(blk(2,:)-1);
                ph = P.phase(ind);
                coeffs = X(ind).*conj(phases(ph+1)); % multiply by conjugate phases
                [s, e] = replab.numerical.sum2(coeffs);
                e = abs(e);
                if any(ph ~= 0 & ph ~= po/2 & ph ~= po/4 & ph ~= 3*po/4)
                    % add the 1e-15 error on phases when relevant
                    % n+1 because there is an error before computing the mean, and error after multiplication of the mean
                    e = e + 1e-15*abs(s)*(n+1);
                end
                s = s/n;
                e = e/n + abs(s)*eps(1/n); % error on 1/n
                E(ind) = e;
                X1(ind) = s * phases(ph+1); % set the coefficients
            end
            if nargout > 1
                eX1 = norm(E, 'fro');
            end
        end

    end

end
