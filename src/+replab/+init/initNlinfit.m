function res = initNlinfit(verbose)
% Verifies that the nlinfit optimization function is available
%
% This loads the ``optim`` Octave package if necessary
%
% Returns:
%   logical: True if the nlinfit function is available
    switch exist('nlinfit')
      case 0
        if replab.compat.isOctave
            try
                if verbose >= 1
                    disp('Loading optim package for Octave...');
                end
                pkg load optim
            catch
                warning('Please install the optim package for Octave to enable compact groups algorithms.');
            end
        else
            warning('Please install the Matlab statistics toolbox to enable continuous groups algorithms.');
        end
      case 2 % all good, it's a .m file (usually Octave)
      case 6 % all good, it's a .p file (usually Matlab)
      otherwise
        error('Strange result %d for exist(''nlinfit'')', exist('nlinfit'));
    end
end
