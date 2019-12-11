classdef ClassElement < replab.Str
    
    properties
        name % charstring: Method or property identifier
        attributes % struct: Attributes from the ``methods``/``properties`` block
        docLines % row cell vector of charstring: Documentation comment lines
    end
    
end
