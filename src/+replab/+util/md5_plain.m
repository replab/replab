function digest = md5_plain(bytes)
% Computes the MD5 hash of a vector of bytes
%
% This follows the pseudo-code on Wikipedia `<https://en.wikipedia.org/wiki/MD5>`_
%
% See `<https://rosettacode.org/wiki/MD5#MATLAB>`_ for another reference MATLAB/Octave
% implementation.
%
% Coding algorithms using binary representation in MATLAB/Octave is painful because
%
% - addition saturates instead of wrapping,
% - integer division rounds to nearest.
%
% Thus the quite convoluted code below! We encode all 32 bits integers in the uint64 type.
%
% All 8 <-> 32 bits conversions are little endian.

    mask32 = uint64(2^32 - 1); % mask for modulo 2^32 addition

    % s specifies the per-round shift amounts
    s = zeros(1, 64);
    s(1:16) = repmat([7 12 17 22], 1, 4);
    s(17:32) = repmat([5 9 14 20], 1, 4);
    s(33:48) = repmat([4 11 16 23], 1, 4);
    s(49:64) = repmat([6 10 15 21], 1, 4);

    % Use binary integer part of the sines of integers (radians) as constants
    K = uint64(floor(abs(sin(1:64)) .* 2^32));

    % Initialize variables
    a0 = uint64(hex2dec('67452301'));
    b0 = uint64(hex2dec('EFCDAB89'));
    c0 = uint64(hex2dec('98BADCFE'));
    d0 = uint64(hex2dec('10325476'));

    bytes = uint64(bytes);
    nBytes = numel(bytes);
    % Pad with 0x00 bytes so that the message length in bytes = 56 (mod 64)
    bytes = [bytes, 128, zeros(1, mod(55 - nBytes, 64))];

    % Convert the message to 32-bit integers
    words = reshape(bytes, 4, numel(bytes) / 4);
    words = words(1,:) + words(2,:)*2^8 + words(3,:)*2^16 + words(4,:)*2^24;

    % Append the bit length as a 64-bit integer (two 32-bits words)
    bitlen = nBytes * 8;
    words = [words, mod(bitlen, 2^32), mod(bitlen / 2^32, 2^32)];
    words = uint64(words); % to be sure

    % Process the message in successive 512-bit chunks
    for j = 1:16:numel(words)
        % Initialize hash value for this chunk
        a = a0;
        b = b0;
        c = c0;
        d = d0;
        % Main loop
        for i = 1:64
            if i <= 16
                % Round 1
                f = bitor(bitand(b, c), bitand(bitcmp(b), d));
                g = i - 1;
            elseif i <= 32
                % Round 2
                f = bitor(bitand(b, d), bitand(c, bitcmp(d)));
                g = mod(5 * i - 4, 16);
            elseif i <= 48
                % Round 3
                f = bitxor(b, bitxor(c, d));
                g = mod(3 * i + 2, 16);
            else
                % Round 4
                f = bitxor(c, bitor(b, uint64(bitcmp(uint32(d)))));
                g = mod(7 * i - 7, 16);
            end
            % Sum modulo 2^32
            f = bitand(a + f + words(j + g) + K(i), mask32);
            % Update a b c d
            tmp = d;
            d = c;
            c = b;
            f = replab.util.bitrol32(f, s(i));
            b = bitand(b + f, mask32);
            a = tmp;
        end
        % Add this chunk's hash to result so far, modulo 2^32
        a0 = bitand(a0 + a, mask32);
        b0 = bitand(b0 + b, mask32);
        c0 = bitand(c0 + c, mask32);
        d0 = bitand(d0 + d, mask32);
    end

    digest = [a0 b0 c0 d0];

    % Convert hash to bytes
    digest = [digest
              bitand(digest, (2^8-1)*2^8) / 2^8
              bitand(digest, (2^8-1)*2^16) / 2^16
              bitand(digest, (2^8-1)*2^24) / 2^24];
    digest = bitand(digest, 2^8-1);
    digest = reshape(digest, 1, numel(digest));

    % Convert bytes to hexadecimal
    digest = dec2hex(digest);
    digest = reshape(transpose(digest), 1, numel(digest));
end
