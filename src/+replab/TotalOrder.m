classdef TotalOrder < replab.Domain
% Describes a domain equipped with a total order

    methods % Comparison

        function c = compare(self, x, y)
        % Compares two elements of a domain
        %
        % The method returns:
        %
        % - ``-1`` if ``x < y``
        % - ``0`` if ``x == y``
        % - ``1`` if ``x > y``
        %
        % Args:
        %   x (element): First element
        %   y (element): Second element
        %
        % Returns:
        %   -1,0,1: Result of comparison
            error('Abstract');
        end

    end

    methods (Access = protected)

        function I = quicksort(self, array, I, lo, hi)
            if lo >= 1 && hi >= 1 && lo < hi
                [I, p] = self.partition(array, I, lo, hi);
                I = self.quicksort(array, I, lo, p - 1);
                I = self.quicksort(array, I, p + 1, hi);
            end
        end

        function [I, i] = partition(self, array, I, lo, hi)
            pivot = array{I(hi)};
            i = lo - 1;
            for j = lo:hi
                if self.compare(array{I(j)}, pivot) <= 0
                    i = i + 1;
                    I([i j]) = I([j i]);
                end
            end
        end

    end

    methods % Sort algorithm

        function I = sort(self, array)
        % Returns a permutation that sorts the given array
        %
        % Returns a permutation ``I`` such that ``array(I)`` is sorted.
        %
        % Args:
        %   array (cell(1,\*) of elements of this domain): Array to sort
        %
        % Returns:
        %   integer(1,\*): Permutation that sorts the array
            I = 1:length(array);
            I = self.quicksort(array, I, 1, length(array));
        end

    end

end
