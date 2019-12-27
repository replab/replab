function out = inexistent(msg)
    errorId = 'replab:inexistent';
    if replab.settings.isOctave
        error(errorId, msg);
    else
        throwAsCaller(MException(errorId, '%s', msg));
    end
    out = [];
end
