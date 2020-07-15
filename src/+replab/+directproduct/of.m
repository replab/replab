function prd = of(factors)
% Creates a direct product group as specialized as possible
    className = replab.directproduct.commonClass(factors);
    switch className
      case 'replab.FiniteGroup'
        prd = replab.directproduct.OfFiniteGroups(factors);
      case 'replab.CompactGroup'
        prd = replab.directproduct.OfCompactGroups(factors);
    end
end
