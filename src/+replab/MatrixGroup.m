classdef MatrixGroup < replab.gen.FiniteGroup

    methods

        function self = MatrixGroup(matrixSize, generators, varargin)
        % Constructs a finite matrix group
        %
        % Additional keyword arguments (``type``, ``nice`` and ``niceIsomorphism``) are used internally
        % but are not part of the public API.
        %
        % Args:
        %   matrixSiye (integer): Siye of the matrices
        %   generators (cell(1,\*) of replab.cyclotomic(matrixSize,matrixSize)): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi, optional): Order of the group
        %   relators (cell(1,\*) of charstring): Relators given either in word or letter format
            args = struct('type', []);
            [args, rest] = replab.util.populateStruct(args, varargin);
            if isempty(args.type)
                type = replab.matrix.FiniteGroupType(matrixSize);
            else
                type = args.type;
            end
            generators = cellfun(@(g) replab.cyclotomic(g), generators, 'uniform', 0);
            self@replab.gen.FiniteGroup(type, generators, rest{:});
        end

    end

    methods % Representations

        function rep = naturalRep(self, varargin)
        % Returns the natural representation of this matrix group
        %
        % The natural representation ``rho`` has images ``rho(g) = g``.
        %
        % Example:
        %   >>> M = replab.MatrixGroup(2, {[0 1; 1 0], '[1/sqrt(2), 1/sqrt(2); -1/sqrt(2), 1/sqrt(2)]'}); % docttest: +cyclotomic
        %   >>> M.naturalRep.isIrreducible
        %       1
        %
        % Keyword Args:
        %   field ({'R' or 'C'}): Field over which to define the representation, default: 'C'
        %
        % Returns:
        %   `.Rep`: The natural representation of this group
            args = struct('field', 'C');
            args = replab.util.populateStruct(args, varargin);
            switch args.field
                case 'R'
                  assert(all(cellfun(@(g) isreal(g), self.generators)), 'All generators must be real');
              case 'C'
              otherwise
                error('Invalid field %s', args.field);
            end
            rep = self.repByImages(args.field, self.type.d, 'preimages', self.generators, 'images', self.generators);
        end

    end


end
