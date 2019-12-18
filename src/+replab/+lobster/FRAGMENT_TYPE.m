classdef FRAGMENT_TYPE
    
% Changed the enumeration to properties as the RepLAB code parsing does not support enumerations well
   properties (Constant)
       VAR = 1
       BLOCK_START = 2
       BLOCK_END = 3
       TEXT = 4
   end
   
end
