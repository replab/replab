function prd = of(H, A)
% Constructs a wreath product with the most refined type
%
% See `+replab.CompactGroup.wreathProduct`
    className = replab.directproduct.commonClass({H A});
    switch className
      case 'replab.NiceFiniteGroup'
        prd = replab.wreathproduct.OfNiceFiniteGroup(H, A);
      case 'replab.FiniteGroup'
        prd = replab.wreathproduct.OfFiniteGroup(H, A);
      case 'replab.CompactGroup'
        prd = replab.wreathproduct.OfCompactGroup(H, A);
    end
end
