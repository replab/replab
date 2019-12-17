%% Pretty printing infrastructure
%
% Matlab has little support for consistent pretty printing of objects. Thus,
% RepLAB contains helper classes and functions to print arbitrary objects in
% a clear style.
%
% Three styles of printing are used by default:
%
% - a 'tiny' style that prints only the size and type of the object,
%
% - a 'short' style that fits on a single display line,
%
% - a 'long' style that can use multiple lines.
%
% The style used when displaying objects in the REPL/command line is the
% 'long' style. For example:
P = replab.Permutations(3)

%%
% which is the output of the 'replab.longStr' function
replab.longStr(P)

%%
% while the short version is much less informative.
replab.shortStr(P)

%%
% All classes in RepLAB inherit the 'replab.Str' base class, which provides
% explicit 'longStr' and 'shortStr' methods. Those methods take a dimension
% limit for width (and possibly height).
% In contrast, the 'replab.shortStr' and 'replab.longStr' use default values
% when those arguments are not given (maxRows = 25, maxColumns = 80).
P.longStr(80, 25)
P.shortStr(80)
%%
% while implementing disp() using 'longStr'.
%
% The default implementation of 'longStr' is to print a short description of
% the object on the first line, followed by public properties.
%
% However, additional fields can be printed by overloading the
% 'additionalFields' method, while properties can be hidden from view by
% overriding the 'hiddenFields' method.
%
% For example, 'FiniteGroup.order' is not a property, as it is computed on
% demand. Thus, when printing a permutation group, the output depends on
% whether the order is known
P = replab.Permutations(3).subgroup({[2 3 1]})
P.order
P
