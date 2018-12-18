function B = isnumeric(INT)
% vpi/isnumeric: true for a vpi array
% usage: B = isnumeric(INT)
% 
% arguments: (input)
%  INT - vpi object or scalar integer numeric values
%
% arguments: (output)
%  B  - scalar true value. true for all vpi input
% 
% Example:
%  INT = vpi([inf 1 -inf NaN]);
%  ans =
%     1
%
%  See also: isnan, isinf
%  
% 
%  Author: John D'Errico
%  e-mail: woodchips@rochester.rr.com
%  Release: 1.0
%  Release date: 2/16/2011


if nargin~=1
  error('VPI:ISNUMERIC:monadicoperator','isnumeric is a monadic function, exactly 1 argument is required')
end

% if we got into here, then by definition, the variable is numeric
B = true;

