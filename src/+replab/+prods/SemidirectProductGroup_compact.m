classdef SemidirectProductGroup_compact < replab.SemidirectProductGroup & replab.CompactGroup

    methods

        function self = SemidirectProductGroup_compact(phi)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            assert(isa(H, 'replab.CompactGroup'));
            assert(isa(N, 'replab.CompactGroup'));
            self.phi = phi;
            self.H = H;
            self.N = N;
            self.identity = {H.identity N.identity};
        end

    end

    methods % Implementations

        % CompactGroup

        function b = hasReconstruction(self)
            b = self.H.hasReconstruction && self.N.hasReconstruction && self.H.reconstruction.source.n == 0;
            % we only support tori in the group acted upon
        end

        function [mu, R] = reconstruction(self)
            [muH, RH] = self.H.reconstruction;
            [muN, RN] = self.N.reconstruction;
            assert(muH.source.n == 0);
            idH = self.H.identity;
            idN = self.N.identity;
            mu = muN.andThen(self.Nembedding);
            setsH = cellfun(@(S) cellfun(@(s) {s, idN}, S, 'uniform', 0), RH.sets, 'uniform', 0);
            setsN = cellfun(@(S) cellfun(@(s) {idH, s}, S, 'uniform', 0), RN.sets, 'uniform', 0);
            R = replab.SetProduct(self, horzcat(setsH, setsN), true);
        end

    end

end
