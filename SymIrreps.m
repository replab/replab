classdef SymIrreps < replab.Str
% Singleton pattern in OOP

    properties
        reps % (cell(1,\*) of whatever): Computed value or [] if not computed
    end

    methods

        function r = rep(self, i)
            if i > length(self.reps) || isempty(self.reps{i})
                self.reps{i} = 1;
                % do the computation
            end
            r = self.reps{i};
        end

    end


    % The part below is about keeping only one instance of SymIrreps
    % in memory; the user will call SymIrreps.instance to grab it
    methods (Access = protected)

        function self = SymIrreps
        % construct
            self.reps = {};
        end

    end

    methods (Static)

        function res = instance
            persistent instance_ % will preserve the value across function calls
                                 % this contains the only copy of SymIrreps
                                 % in memory
            if isempty(instance_)
                res = SymIrreps; % call the constructor
            end
        end

    end

end
