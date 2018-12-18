classdef Domain < replab.cat.Laws
    
    properties (Abstract, SetAccess = protected)
        canEqv;
        canHash;
        canSample;
    end
    
    methods
        
        % Implement
        %
        % function b = eqv(self, x, y)
        % function h = hash(self, x)
        % function x = sample(self)
        
    end
    
    methods
        
        function law_eqv_D(self, x)
            self.assertTrue(self.eqv(x, x));
        end
        
        function law_hash_eqv_DD(self, x, y)
            if self.hash(x) ~= self.hash(y)
                self.assertNotEqv(x, y);
            end
        end
        
        function assertNotEqv(self, x, y, context)
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
            
            message = replab.cat.lawsMessage(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end            

        end
            
        function assertEqv(self, x, y, context)
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
            
            message = replab.cat.lawsMessage(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end            

        end

    end
    
    methods (Static)
        
        function h = hashIntegers(M)
            if length(M) == 1
                h = double(31 + int32(M));
            else
                h = javaMethod('hashCode', 'java.util.Arrays', int32(M(:)));
            end
        end
        
        function D = integerRange(minValue, maxValue)
            desc = sprintf('Integers in [%d, ..., %d]', minValue, maxValue);
            sampleFun = @() randi([minValue maxValue]);
            hashFun = @(x) x;
            D = replab.cat.DomainFun(desc, @isequal, hashFun, sampleFun);
        end
        
        function D = signedIntegers(n)
            desc = sprintf('Integers {-%d, ..., -1, 1, ..., %d}', n, n);
            sampleFun = @() randi(n) * (randi([0 1])*2-1);
            hashFun = @(x) x;
            D = replab.cat.DomainFun(desc, @isequal, hashFun, sampleFun);
        end

        function D = integerColumnVectors(d, minValue, maxValue)
            desc = sprintf('Column vectors of length %d with coefficients in [%d, ..., %d]', minValue, maxValue);
            sampleFun = @() randi([minValue, maxValue], d, 1);
            hashFun = @(M) replab.cat.Domain.hashIntegers(M);
            D = replab.cat.DomainFun(desc, @isequal, hashFun, sampleFun);
        end
        
        function D = integerMatrices(m, n, minValue, maxValue)
            desc = sprintf('Integer %d x %d matrices with coefficients in [%d, ..., %d]', m, n, minValue, maxValue);
            sampleFun = @() randi([minValue maxValue], m, n);
            hashFun = @(M) replab.cat.Domain.hashIntegers(M);
            D = replab.cat.DomainFun(desc, @isequal, hashFun, sampleFun);
        end
        
    end
    
end
