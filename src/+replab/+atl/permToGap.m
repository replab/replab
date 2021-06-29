        function s = permToGap(g)
        % Write a GAP System string representation of a RepLAB permutation
        %
        % Args:
        %   g (permutation): Permutation
        %
        % Returns:
        %   charstring: String representation
            s = sprintf('Inverse(PermList([%s]))', strjoin(arrayfun(@num2str, g, 'uniform', 0), ','));
        end
