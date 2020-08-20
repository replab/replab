classdef Group < replab.Monoid
% Describes a group
%
% A group is a `.Monoid` where each element has an inverse.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.GroupLaws(self);
        end

    end

    methods % Group methods

        % Abstract methods

        function xInv = inverse(self, x)
        % Computes the inverse of an element
        %
        % Given ``x``, returns ``xInv`` such that
        %
        % ``x xInv = identity``
        %
        % Args:
        %   x (element): Group element to compute the inverse of
        %
        % Returns:
        %   element: Inverse of ``x``
            error('Abstract');
        end

        % Methods with default implementations

        function x = leftConjugate(self, by, on)
        % Returns the left conjugate of a group element
        %
        % Convenience method that can be overriden for speed optimizations
        %
        % Args:
        %   by (element): Element conjugating
        %   on (element): Element conjugated
        %
        % Returns:
        %   element: left conjugate, i.e. ``by * on * by^-1`` in multiplicative notation
            x = self.composeWithInverse(self.compose(by, on), by);
        end

        function z = commutator(self, x, y, startWithInverse)
        % Returns the commutator of two group elements
        %
        % Args:
        %   x (element): First element
        %   y (element): Second element
        %   startWithInverse (logical): True if we compute ``x^-1 y^-1 x y``, false if we compute ``x y x^-1 y^-1``
        %
        % Returns:
        %   element: Commutator
            xy = self.compose(x, y);
            yx = self.compose(y, x);
            if startWithInverse
                xIyI = self.inverse(yx);
                z = self.compose(xIyI, xy);
            else
                z = self.composeWithInverse(xy, yx);
            end
        end

        function z = composeWithInverse(self, x, y)
        % Returns the composition of an element with the inverse of another element
        %
        % Convenience method that can be overriden for speed optimizations
        %
        % Args:
        %   x (element): First element
        %   y (element): Second element
        %
        % Returns:
        %   element: the result of ``x * y^-1`` in multiplicative notation
            z = self.compose(x, self.inverse(y));
        end

        function z = composeLetters(self, array, letters)
        % Returns the composition of elements from an array picked according to the given indices/letters
        %
        % Computes ``composeAll({array(letters(1)) ... array(letters(n))})`` with the convention that
        %
        % Args:
        %   letters (integer(1,\*)): Letters composing word
        %   array (cell(1,\*) of elements): Elements to compose
        %
        % Returns:
        %   element: The result of the composition
            if nargin < 4
                arrayInverse = [];
            end
            z = self.identity;
            L = length(letters);
            for i = 1:L
                l = letters(i);
                if l > 0
                    z = self.compose(z, array{l});
                else
                    z = self.composeWithInverse(z, array{-l});
                end
            end
        end

    end

    methods % Morphisms

        function G = automorphismGroup(self)
        % Returns the automorphism group of this group
            G = replab.AutomorphismGroup(self);
        end

        function m = innerAutomorphism(self, by)
        % Returns the inner automorphism given using conjugation by the given element
        %
        % Args:
        %   by (element): Conjugating element
        %
        % Returns:
        %   `.Morphism`: Group automorphism
            byInv = self.inverse(by);
            m = replab.Isomorphism.lambda(self, self, @(h) self.leftConjugate(byInv, h), @(h) self.leftConjugate(by, h));
        end

    end

    methods (Static)

        function group = lambda(header, eqvFun, sampleFun, composeFun, identity, inverseFun)
        % Constructs a group from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the parent `.Domain.eqv` method
        %   sampleFun (function_handle): Handle implementing the parent `.Domain.sample` method
        %   composeFun (function_handle): Handle implementing the parent `.Monoid.compose` method
        %   identity (element): Identity element of this monoid
        %   inverseFun (function_handle): Handle implementing the `inverse` method
        %
        % Returns:
        %   `+replab.Group`: The constructed group

            group = replab.lambda.Group(header, eqvFun, sampleFun, ...
                                        composeFun, identity, inverseFun);
        end

    end

end
