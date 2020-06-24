function string = seconds2human(secs, outputFormat)
% Converts the given number of seconds into a human-readable string.
%
% Returns a human-readable string from a given (usually large) amount of seconds.
%
% Example:
%   >>> str = replab.infra.repl.seconds2human(1463456.3)
%       str =
%       'About 2 weeks and 2 days'
%
%
% You may also call the function with a second input argument; either
% ``'short'`` (the default) or ``'full'``. This determines the level of detail
% returned in the string:
%
% Example:
%   >>> str = replab.infra.repl.seconds2human(1463456.3, 'full')
%       str =
%       '2 weeks, 2 days, 22 hours, 30 minutes, 56 seconds'
%
%
% The ``'short'`` format returns only the two largest units of time.
%
% NOTE: `seconds2human` defines one month as an "average" month, which
% means that the string 'month' indicates 30.471 days.
%
% Args:
%   secs (double): Time in seconds
%   outputFormat({'short', 'full'}, optional): Short or full output format, default ``'short'``
%
%
% Returns:
%   charstring: Formatted duration
%
% Note:
%   Adapted for RepLAB use from the source code by:
%
%   Name   : Rody P.S. Oldenhuis
%   E-mail : oldenhuis@gmail.com
%   Licence: 2-clause BSD (See Licence.txt)
%   If you find this work useful, please consider a donation:
%   https://www.paypal.me/RodyO/3.5
%
%   We changed the string ending from the period ``'.'`` to nothing,
%   and removed the code that handles matrices.

    strEnding = ''; % replace by '.' to get old behavior

    assert(isscalar(secs) && isa(secs, 'double'));

    % define some intuitive variables
    Seconds   = round(1                 );
    Minutes   = round(60     * Seconds  );
    Hours     = round(60     * Minutes  );
    Days      = round(24     * Hours    );
    Weeks     = round(7      * Days     );
    Months    = round(30.471 * Days     );
    Years     = round(365.26 * Days     );
    Centuries = round(100    * Years    );
    Millennia = round(10     * Centuries);

    % put these into an array, and define associated strings
    units   = [Millennia, Centuries, Years, Months, Weeks, ...
               Days, Hours, Minutes, Seconds];
    singles = {'millennium'; 'century'; 'year'; 'month'; ...
               'week'; 'day'; 'hour'; 'minute'; 'second'};
    plurals = {'millennia' ; 'centuries'; 'years'; 'months'; ...
               'weeks'; 'days'; 'hours'; 'minutes'; 'seconds'};

    % cut off all decimals from the given number of seconds
    assert(isnumeric(secs), 'seconds2human:seconds_mustbe_numeric', ...
           'The argument ''secs'' must be a scalar or matrix.');
    secs = round(secs);

    % parse second argument
    if nargin < 2
        outputFormat = 'short';
    end
    assert(ischar(outputFormat), 'The second argument must be either ''short'' or ''full''.');
    switch outputFormat
      case 'full'
        short = false;
      case 'short'
        short = true;
      otherwise
        error('The second argument must be either ''short'' or ''full''.');
    end

    counter = 0;
    if short
        string = 'About ';
    else
        string = '';
    end

    % possibly quick exit
    if (secs < 1)
        string = ['Less than one second' strEnding];
        return
    end

    % build string for j-th amount of seconds
    for i = 1:length(units)
        % amount of this unit
        amount = fix(secs/units(i));
        % include this unit in the output string
        if amount > 0
            % increase counter
            counter = counter + 1;
            % append (single or plural) unit of time to string
            if (amount > 1)
                string = [string, num2str(amount), ' ', plurals{i}];%#ok
            else
                string = [string, num2str(amount), ' ', singles{i}];%#ok
            end
            % Finish the string after two units if short format is requested
            if (counter > 1 && short)
                string = [string strEnding];
                return
            end %#ok

            % determine whether the ending should be strEnding or a comma (,)
            if (rem(secs, units(i)) > 0)
                if short
                    ending = ' and ';
                else
                    ending = ', ';
                end
            else
                ending = strEnding;
            end
            string = [string, ending]; %#ok
        end
        % subtract this step from given amount of seconds
        secs = secs - amount*units(i);
    end

end % seconds2human
