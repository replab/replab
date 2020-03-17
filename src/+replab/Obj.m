classdef Obj < replab.Str
% Base class that provides sane pretty printing and instance equality tests

    methods (Static)

        function id = newId
            persistent lastId
            if isempty(lastId)
                lastId = uint64(0);
            end
            id = lastId + uint64(1);
            lastId = id;
        end

    end

    properties (Access = protected)
        id_ % (uint64): Unique object ID
    end

    methods


        function i = id(self)
            if isempty(self.id_)
                self.id_ = replab.Obj.newId;
            end
            i = self.id_;
        end

        function res = eq(self, rhs)
        % Equality test
        %
        % Workaround bug of == not implemented for handles
        %
        % Args:
        %   self (object): first object
        %   rhs (object): second object to compare to
        %
        % Returns:
        %   boolean: true iff self == rhs
            if ~isa(rhs, 'replab.Obj')
                res = false;
                return
            end
            l = arrayfun(@(x) x.id, self);
            r = arrayfun(@(x) x.id, rhs);
            res = (l == r);
        end

    end

end
