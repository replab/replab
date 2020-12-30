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

end
