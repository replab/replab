function out = inexistent(msg)
    errorId = 'replab:inexistent';
    if replab.platformIsOctave()
        error(errorId, msg);
    else
        throwAsCaller(MException(errorId, '%s', msg));
    end
    out = [];
end
