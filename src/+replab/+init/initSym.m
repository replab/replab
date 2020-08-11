function res = initSym
% Verifies that symbolic computation is available
%
% This loads the ``symbolic`` Octave package if necessary
%
% Returns:
%   logical: True if the Symbolic Toolbox is available

    verbose = replab.globals.verboseInit;
    res = false;
    switch exist('syms')
      case 0
        if replab.compat.isOctave
            try
                if verbose >= 1
                    disp('Loading symbolic package for Octave...');
                end
                pkg load symbolic
            catch
                warning('Please install the symbolic package to enable exact .repByImages');
            end
        else
            warning('Please install the Matlab Symbolic toolbox to enable exact .repByImages');
        end
      case 2 % all good, it's a .m file (recent Octave and Matlab)
        res = true;
      case 6 % all good, it's a .p file (not seen)
        res = true;
      otherwise
        error('Strange result %d for exist(''syms'')', exist('syms'));
    end
end
