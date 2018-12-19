classdef SCProjections < replab.cat.Laws
    
    methods
        
        
        function y = sampleGaussianWithVariances(self, weights)
        % weight is a 1 x N vector of integers
        % such y(i) is normally distributed with mean 0
        % and variance 2/weight(i)
            error('Not implemented');
        end
       
        function y = sampleSelfAdjointCoefficients(self, sc)
            error('Not implemented');
        end
        
        function y = sampleCoefficients(self, sc)
            error('Not implemented');
        end
        
        function y = projectCoefficients(self, x)
        end
        
        function y = projectSelfAdjointCoefficients(self, sc, x)
        end
        
        function x = toSelfAdjointMatrix(y)
        end
        
        function x = toMatrix(y)
        end
        
    end
    
end
