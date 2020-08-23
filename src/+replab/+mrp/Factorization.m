classdef Factorization < replab.Obj
% Computes the factorization of a permutation group elements in its generators

    properties (SetAccess = protected)
        group % (`+replab.PermutationGroup`): Group which elements are factorized
    end

    methods

        function letters = factorize(self, g)
        % Returns the letters that compose the word with the given image
        %
        % Args:
        %   g (permutation in `.group`): Permutation to factorize
        %
        % Returns:
        %   integer(1,\*): Corresponding word expressed in letters
            error('Abstract');
        end

    end

end
