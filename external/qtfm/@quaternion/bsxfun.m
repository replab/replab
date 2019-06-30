function C = bsxfun(fun, A, B)
% BSXFUN  Binary Singleton Expansion Function
% (Quaternion overloading of standard Matlab function.)

% The reason for providing this function (apart from its possible intrinsic
% usefulness) is that many Matlab functions such as var, cov, std, kron,
% depend on bsxfun. Early releases of QTFM worked with these functions at
% the time, but as Matlab was changed to use bsxfun they failed to work
% with QTFM, because no overloading of bsxfun was provided.

% Copyright (c) 2009, 2010, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% Rewritten completely in 2015.
% See the file : Copyright.m for further details.

narginchk(3, 3), nargoutchk(0, 1)

if ~isa(fun, 'function_handle')
    error('First parameter must be a function handle')
end

% We do not check that the function handle supplied is meaningful, since
% the algorithm implemented here as of October 2015 will work for any
% meaningful function that obeys the rules for bsxfun (see the
% documentation on the Matlab version of bsxfun for these rules).

% We don't check the types of the A and B parameters either. They could
% both be quaternions, or one could be numeric/logical. The function
% identified by the handle fun will fail with an error if it can't handle
% the mix of types, so we leave the error checking to fun, rather than make
% this code overcomplicated. Similarly for quaternions with different
% element types (e.g. double and int8), we leave the checking to fun.

CS = zeros(1, max(ndims(A), ndims(B))); % For the size of the output array.
AS = zeros(size(CS)); % For the size of the A input.
BS = zeros(size(CS)); % For the size of the B input.

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

SA = AS == 1; % Logical array indicating singleton dimensions of A.
SB = BS == 1; % Logical array indicating singleton dimensions of B.

% Pre-allocate the output array. This array could be logical or quaternion,
% depending on the function fun.

if any(strcmp(func2str(fun), {'eq', 'ne', 'gt', 'ge', 'lt', 'le'}))
    C = false(CS);
else
    C = zerosq(CS, class(A.x));
end

% Every element of C must have a matching element in each of A and B. There
% are no singleton dimensions in C. Therefore we can use linear indices
% into C, and for each one compute a set of subscripts according to the
% number of dimensions in C. We can then modify these subscripts to allow
% for singleton dimensions in A or B, then compute linear indices into A
% and B. Once we have those we can apply the function fun to A and B using
% the linear indices and store results into C using linear indices. Clearly
% this won't work in one go if C is large, because the linear index arrays
% will also be large. In this case, we need to process C in chunks of
% linear indices, splitting it into manageable pieces.

% Initial implementation, processing all of C in one go.

IC = (1:prod(CS))'; % Linear indices into C. We carefully make this a
                    % column vector, so that all the subscript arrays are
                    % returned as column vectors. If we don't do this, some
                    % may be row vectors, and we'll get errors.
                    
% To be done in chunks, ultimately.

switch ndims(C)
    case 1
        % This cannot happen, since even a vector has two dimensions (one
        % of which will be unity).

        error('Program error in bsxfun: ndims returned 1.')
    case 2
        [C1, C2] = ind2sub(CS, IC); % Compute subscripts into C.
        
        IA = sub2ind(AS, sex(SA(1), C1), sex(SA(2), C2)); % Indices into A.
        IB = sub2ind(BS, sex(SB(1), C1), sex(SB(2), C2)); % Indices into B.
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
