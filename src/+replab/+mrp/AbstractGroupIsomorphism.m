classdef AbstractGroupIsomorphism < replab.FiniteIsomorphism

    methods

        function iso1 = withUpdatedSource(self, source1)
        % Returns a copy of this isomorphism with the source group modified
        %
        % The new source group and the original source group must only differ in the generator names,
        % otherwise undefined results will occur.
        %
        % Args:
        %   source1 (`+replab.AbstractGroup`): New source
        %
        % Returns:
        %   `.AbstractGroupIsomorphism`: Updated isomorphism
            error('Abstract');
        end

    end

end
