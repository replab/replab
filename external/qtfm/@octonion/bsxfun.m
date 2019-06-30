function C = bsxfun(fun, A, B)
% BSXFUN  Binary Singleton Expansion Function
% (Octonion overloading of standard Matlab function.)

% The reason for providing this function (apart from its possible intrinsic
% usefulness) is that many Matlab functions such as var, cov, std, kron,
% depend on bsxfun. Early releases of QTFM worked with these functions at
% the time, but as Matlab was changed to use bsxfun they failed to work
% with QTFM, because no overloading of bsxfun was provided.

% This code is a copy of the quaternion code, with very minor changes.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 1)

if ~isa(fun, 'function_handle')
    error('First parameter must be a function handle')
end

CS = zeros(1, max(ndims(A), ndims(B)));
AS = zeros(size(CS));
BS = zeros(size(CS));

for D = 1:length(CS)
    
    % Check that the sizes of the two arrays are compatible in dimension D.

    AS(D) = size(A, D);
    BS(D) = size(B, D);
    
    CS(D) = max(AS(D), BS(D));
    
    if AS(D) == BS(D) || AS(D) == 1 || BS(D) == 1
        continue
    else
        % This is the error message that Matlab's bsxfun gives:
        error(['Non-singleton dimensions of the two input arrays ', ...
               'must match each other.'])
    end
end

SA = AS == 1;
SB = BS == 1;

if any(strcmp(func2str(fun), {'eq', 'ne', 'gt', 'ge', 'lt', 'le'}))
    C = false(CS);
else
    C = zeroso(CS, class(part(A, 1)));
end

IC = (1:prod(CS))'; 

switch ndims(C)
    case 1
        % This cannot happen, since even a vector has two dimensions (one
        % of which will be unity).

        error('Program error in bsxfun: ndims returned 1.')
    case 2
        [C1, C2] = ind2sub(CS, IC);
        
        IA = sub2ind(AS, sex(SA(1), C1), sex(SA(2), C2));
        IB = sub2ind(BS, sex(SB(1), C1), sex(SB(2), C2));
    case 3
        [C1, C2, C3] = ind2sub(CS, IC);
        
        IA = sub2ind(AS, sex(SA(1), C1), ...
                         sex(SA(2), C2), ...
                         sex(SA(3), C3));
        IB = sub2ind(BS, sex(SB(1), C1), ...
                         sex(SB(2), C2), ...
                         sex(SB(3), C3));
    case 4
        [C1, C2, C3, C4] = ind2sub(CS, IC);
        
        IA = sub2ind(AS, sex(SA(1), C1), ...
                         sex(SA(2), C2), ...
                         sex(SA(3), C3), ...
                         sex(SA(4), C4));
        IB = sub2ind(BS, sex(SB(1), C1), ...
                         sex(SB(2), C2), ...
                         sex(SB(3), C3), ...
                         sex(SB(4), C4));

    otherwise
        error('bsxfun cannot yet handle more than four dimensions')
end

% Since this is a class method, normal array indexing does not work.
% See the file implementation_notes.txt, point 8.

% C(IC) = fun(A(IA), B(IB)); % This is what we are doing below.

C = subsasgn(C, substruct('()', {IC}), ...
    fun(subsref(reshape(A, prod(AS), 1), substruct('()', {IA})), ...
        subsref(reshape(B, prod(BS), 1), substruct('()', {IB}))));
end

function O = sex(L, I)
% Singleton expansion (in effect). Given a logical value L, returns an
% array based on I, containing either the same values as I (L == false) or
% ones (L == true). TODO Consider how to eliminate the actual expansion
% here and replace it with 'virtual' expansion.
if L
    O = ones(size(I));
else
    O = I;
end
end

% $Id: bsxfun.m 1004 2017-11-15 17:14:09Z sangwine $
