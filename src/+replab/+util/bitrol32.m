function res = bitrol32(num, shift)
% Rotates to the left the given binary number
%
% The binary representation is split into ``[left part][right part]`` where
% - left part is of length "shift"
% - right part is of length "32-shift"
%
% Then the result is simply ``[right part][left part]``
%
% Args:
%   num (uint32 or uint64): Number to shift
%   shift (0..31): Number of positions to shift
    leftMask = (2^shift-1)*2^(32-shift);
    leftQuotient = 2^(32-shift);
    leftRes = bitand(num, leftMask)/leftQuotient;
    rightMask = 2^(32-shift)-1;
    rightFactor = 2^shift;
    rightRes = bitand(num, rightMask)*rightFactor;
    res = leftRes + rightRes;
end
