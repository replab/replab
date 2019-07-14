classdef PermutationBSGSGroup < replab.PermutationGroup & replab.BSGSGroup
    
    methods
        
        function self = PermutationBSGSGroup(parent, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            self@replab.PermutationGroup(parent.domainSize);
            self.action = parent.naturalAction;
            self.generators = generators;
            self.order_ = orderOpt;
        end
        
        function g = sample(self)
        % Repeating method in BSGSGroup to workaround an Octave
            if self.knownChain
                g = self.sampleUniformly;
            else
                g = self.randomBag.sample;
            end
        end
        
    end
    
end
