classdef ConjugacyClass < replab.ConjugacyClass & replab.gen.FiniteSet

    methods

        function self = ConjugacyClass(type, nice, niceIsomorphism, varargin)
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
            args = struct('group', []);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.group)
                self.group = niceIsomorphism.preimageGroup(nice.group);
            else
                self.group = args.group;
            end
        end

    end

    methods % Implementations

        % ConjugacyClass

        function o = elementOrder(self)
            o = self.nice.elementOrder;
        end

        function l = knownRepresentativeCentralizer(self)
            l = self.nice.knownRepresentativeCentralizer;
        end

        function G = representativeCentralizer(self)
            G = self.niceIsomorphism.preimageGroup(self.nice.representativeCentralizer);
        end

    end

end
