classdef IndexedFamily < replab.Domain
% Describes an indexed family of elements
%
% See https://en.wikipedia.org/wiki/Indexed_family , and our indices are
% (bounded) integers, represented by ``vpi`` instances.
%
% The family supports element indexing and searching for elements.

    properties (SetAccess = protected)
        nElements % vpi: Number of elements contained in this enumerator
    end

    methods

        function obj = at(self, ind)
        % Retrieves a element by position
        %
        % Normally, RepLAB encodes big integers using the ``vpi`` type. As a user convenience
        % (esp. on the command line), this method must also accept string and double arguments.
        %
        % Args:
        %   ind (vpi or double or string): Index of element to be retrieved, 1 <= ``ind`` <= ``self.nElements``
        %
        % Returns:
        %   The element at the "ind" position
            error('Abstract');
        end

        function ind = find(self, obj)
        % Returns the index of a given element
        %
        % If the element is not part of the family, this returns ``vpi(0)``.
        %
        % Args:
        %   obj: Element to retrieve
        %
        % Returns:
        %   vpi: 1-based index of the given element
            error('Abstract');
        end

        function C = toCell(self)
        % Returns a row cell array containing all elements of this family
        %
        % Returns:
        %   row cell array of elements: A cell array ``C`` such that C{i} = self.at(i)
        %
        % Raises:
        %   An error if the enumerator is too big and the elements cannot fit in a cell array.
            if self.nElements == 0
                C = cell(1, 0);
            else
                n = self.nElements;
                msg = 'Indexed family of size %s too big to enumerate in a matrix';
                assert(n < intmax('int32'), msg, num2str(self.nElements));
                n = double(n);
                C = cell(1, n);
                for i = 1:n
                    C{i} = self.at(i);
                end
            end
        end

    end

    methods % Implementations

        % Str

        function s = shortStr(self, maxColumns)
            s = sprintf('Indexed family of %s elements', replab.shortStr(self.nElements, maxColumns));
        end

        function lines = longStr(self, maxRows, maxColumns)
            if self.nElements > maxRows - 1
                n = maxRows - 2;
                start = ceil(n/2);
                finish = n - start;
                omit = self.nElements - start - finish;
                omitting = 1;
            else
                start = double(self.nElements);
                omit = [];
                omitting = 0;
                finish = 0;
            end
            table = cell(start + omitting + finish, 3);
            ind = 1;
            for i = 1:start
                table{ind,1} = replab.shortStr(i, maxColumns);
                table{ind,2} = ' = ';
                table{ind,3} = replab.shortStr(self.at(i), maxColumns);
                ind = ind + 1;
            end
            if omitting
                table{ind,1} = ['.. ' replab.shortStr(omit, maxColumns)];
                table{ind,2} = '';
                table{ind,3} = 'elements omitted';
                ind = ind + 1;
            end
            for i = 1:finish
                index = self.nElements - finish + i;
                table{ind,1} = replab.shortStr(index, maxColumns);
                table{ind,2} = ' = ';
                table{ind,3} = replab.shortStr(self.at(index), maxColumns);
                ind = ind + 1;
            end
            lines = vertcat({self.shortStr(maxColumns)}, replab.str.align(table, 'rcl'));
        end

        % Obj

        function l = laws(self)
            l = replab.laws.IndexedFamilyLaws(self);
        end

        % Domain

        function res = eqv(self, x, y)
            indX = self.at(x);
            indY = self.at(y);
            res = indX == indY;
        end

        function obj = sample(self)
            obj = self.at(randint(self.nElements)); % use randint as it is the method equivalent to randi on @vpi
        end

    end

    methods (Static) % IndexedFamily construction

        function family = lambda(nElements, atFun, findFun)
        % Constructs an indexed family from function handles
        %
        % Args:
        %   nElements (vpi): Size of the indexed family, so that the index set is 1..``nElements``
        %   atFun (function_handle): Handle that implements the ``at`` method
        %                            To simplify implementation, it is guaranteed
        %                            that ``atFun`` will receive an argument of type ``vpi``.
        %   findFun (function_handle): Handle that implements the ``find`` method
        %
        % Returns:
        %   `.IndexedFamily`: The constructed indexed family
            family = replab.lambda.IndexedFamily(nElements, atFun, findFun);
        end

    end

end
