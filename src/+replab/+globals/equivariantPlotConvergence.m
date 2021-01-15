function flag = equivariantPlotConvergence(newValue)
% Whether to plot the convergence of equivariant samples
    persistent flagValue
    if nargin == 1
        flagValue = newValue;
    end
    if isempty(flagValue)
        flag = false;
    else
        flag = flagValue;
    end
end
