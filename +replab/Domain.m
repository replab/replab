classdef Domain < replab.Str
% Describes a set of elements with a common structure
%
% Those elements can be compared ("eqv"), and random elements
% can be produced ("sample")
    
    methods % ABSTRACT
        
        function b = eqv(self, t, u)
        % Returns true when t and u are equivalent, and false otherwise
            f = self.eqvFun;
            b = f(t, u);
        end
        
        function t = sample(self)
        % Returns a random element sampled from this domain;
        % no particular guarantees are made about genericity, this
        % method is used for law checks
            f = self.sampleFun;
            t = f();
        end
        
    end
    
    methods % Test helpers
       
        function assertNotEqv(self, x, y, context)
        % Asserts that "x" and "y" are equivalent
            if self.eqv(x, y)
                errorDesc = 'The values %s and %s are equivalent, but they should not be';
                errorId = 'assertNotEqual:equal';
            else
                return
            end
            
            if nargin < 4
                context = '';
            end
            
            names = evalin('caller', 'who');
            nV = length(names);
            values = cell(1, nV);
            for i = 1:nV
                values{i} = evalin('caller', names{i});
            end
            
            message = replab.laws.message(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end            

        end
            
        function assertEqv(self, x, y, context)
        % Asserts that "x" and "y" are not equivalent
            if ~self.eqv(x, y)
                errorDesc = 'The values %s and %s are not equivalent';
                errorId = 'assertEqual:nonEqual';
            else
                return
            end
            
            if nargin < 4
                context = '';
            end
            
            names = evalin('caller', 'who');
            nV = length(names);
            values = cell(1, nV);
            for i = 1:nV
                values{i} = evalin('caller', names{i});
            end
            
            message = replab.laws.message(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end            

        end

    end
    
end
