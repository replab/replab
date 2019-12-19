classdef OnDemandSamples < replab.Str

    properties (SetAccess = protected)
        rep % replab.Rep: Representation to compute things from
        commutantSamples % row cell array of double matrix: Samples from the commutant
        trivialSamples % row cell array of double matrix: Self adjoint samples (TODO ?) from the trivial space
    end

    properties (Access = protected)
        trivial_ % replab.Equivariant: Equivariant space from rep to the trivial representation of size rep.dimension
    end
    
    methods

        function self = OnDemandSamples(rep)
            self.rep = rep;
            self.commutantSamples = {};
            self.trivialSamples = {};
        end

        function T = trivial(self)
            if isempty(self.trivial_)
                trivialRep = self.rep.group.trivialRep(self.rep.field, self.rep.dimension);
                self.trivial_ = trivialRep.equivariant(self.rep);
            end
            T = self.trivial_;
        end
        
        function X = commutantSample(self, i)
            while length(self.commutantSamples) < i
                self.commutantSamples{1, end+1} = self.rep.commutant.sample;
            end
            X = self.commutantSamples{i};
        end

        function X = trivialSample(self, i)
            while length(self.trivialSamples) < i
                self.trivialSamples{1, end+1} = self.trivial.sample;
            end
            X = self.trivialSamples{i};
        end

    end

end
