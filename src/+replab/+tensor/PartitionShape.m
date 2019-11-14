classdef PartitionShape
    
    properties (SetAccess = immutable)
        partition % `+replab.Partition`: Blocks of indices
        shapes % row cell array of `+replab.+tensor.Shape`: Shapes corresponding to the partition blocks
    end
    
    methods
       
        function self = PartitionShape(partition, shapes)
            
        end
        
    end

end
