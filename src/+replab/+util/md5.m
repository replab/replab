function res = md5(bytes)
% Returns the MD5 hash of an ASCII string
%
% Example:
%   >>> replab.util.md5('The quick brown fox jumps over the lazy dog')
%      '9E107D9D372BB6826BD81D3542A419D6'
%   >>> replab.util.md5('')
%      'D41D8CD98F00B204E9800998ECF8427E'
%
% Args:
%   bytes (charstring): Bytes to hash
%
% Returns:
%   charstring: MD5 hash in hexadecimal notation
    if usejava('jvm')
        md = javaMethod('getInstance', 'java.security.MessageDigest', 'MD5');
        res = typecast(md.digest(uint8(bytes)), 'uint8');
        res = strjoin(arrayfun(@(b) dec2hex(b, 2), res, 'uniform', 0), '');
    elseif replab.compat.isOctave
        res = upper(hash('md5', bytes));
    else
        res = replab.util.md5_plain(bytes);
    end
end
