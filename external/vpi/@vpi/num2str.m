function str = num2str(N,varargin)
% vpi/num2str: converts VPI integers to character string form
% usage: str = num2str(N)
% 
% arguments: (input)
%  N - vpi number
%
% arguments: (output)
%  str - character string form for N

if nargin > 1
  warning('vpi:NUM2STR:TBD','Two argument num2str is not yet implemented.')
end

str = disp(N);

