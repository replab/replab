classdef Trilean < replab.Domain
% Describes three-valued logic
%
% See https://en.wikipedia.org/wiki/Three-valued_logic
%
% True and false are represented by `logical` values, while unknown
% is represented by ``[]``.
    
    methods
        
        %% Str methods
        
        function s = headerStr(self)
            s = 'Trilean logic values';
        end
        
        %% Domain methods
        
        function b = eqv(self, t1, t2)
            b = isequal(t1, t2);
        end
        
        function t = sample(self)
            switch randi(3)
              case 1
                t = true;
              case 2
                t = false;
              case 3
                t = [];
            end
        end
                
    end

    methods (Static)
       
        function res = and(values)
        % Returns the 'and' of trilean values
        %
        % Args:
        %   values (row cell array): Trilean values
        %
        % Returns:
        %   true or false or []: Result of the 'and' operation
            if isempty(values)
                res = true; % neutral element for and
                return
            end
            isUnknown = cellfun(@(x) isempty(x), values);
            if any(isUnknown)
                res = [];
                return
            end
            res = all(cell2mat(values));
        end
        
        function res = or(values)
        % Returns the 'or' of trilean values
        %
        % Args:
        %   values (row cell array): Trilean values
        %
        % Returns:
        %   true or false or []: Result of the 'or' operation
            if isempty(values)
                res = false; % neutral element for or
                return
            end
            isUnknown = cellfun(@(x) isempty(x), values);
            if any(isUnknown)
                res = [];
                return
            end
            res = any(cell2mat(values));
        end
        
    end
end
