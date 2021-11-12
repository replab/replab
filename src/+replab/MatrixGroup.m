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
            self@replab.gen.FiniteGroup(type, generators, rest{:});
        end

    end

end
