% Configuration script
%
% Run during the first ``replab_init`` call; or when ``replab_init`` is called after a ``clear all``.
%
% This is essential for big integer suppoort (used everywhere to support finite groups)
replab.init.vpi().require;

% Rich help system; overloads the Matlab/Octave 'help' function. Comment out if incompatible.
% This sets up the `+replab.+globals.defaultHelpFunction` and `+replab.+globals.runsInMatlabEmacs`
% global variables, that are implemented using ``persistent`` variables in the functions. The functions
% are locked with ``mlock`` to avoid being cleared when ``clear all`` is used.
replab.init.initHelp;

% Support for the convex modeling framework and SDP solvers
replab.init.YALMIP().require;
if ~replab.init.existingSdpSolver().works
    replab.init.log_(2, 'Trying to use the embedded solver');
    replab.init.sdpt3().require;
end

% This can be commented out if you're not running tests
replab.init.MOcov().require;
replab.init.MOxUnit().require;

replab.init.cyclolab().require;

replab.init.intlab().require;

% Default values for coset enumeration parameters
replab.globals.cosetEnumerationMethod('R'); % try 'C' if the method doesn't work
replab.globals.maxCosets(2^22);
replab.globals.maxDeductions(100);

% Default values for group recognition
replab.globals.fastChainDomainSize(1000); % only attempt for permutation realizations of domain size <= value
replab.globals.fastChainOrder(10000); % only attempt for groups of order <= value

% Set YOLO mode off and verified arithmetic on
replab.globals.yolo(false);
replab.globals.verifiedArithmetic(true);

% Read additional JSON atlas entries
replab.Atlas.readFolder(fullfile(replab.globals.replabPath, 'atlas'));

% If you have GAP 4 installed, set the path below
% replab.globals.gapBinaryPath('/opt/gap4/bin/gap.sh');
