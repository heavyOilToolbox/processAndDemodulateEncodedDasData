function [Y, fc, fs, INI_PHASE, CARRAMP, NUM, DEN] = ParseInputArguments(Y, fc, fs, varargin)
    if nargin < 3
        error('Not Enough Input Arguments')
    elseif nargin == 3
        INI_PHASE = 0;
        CARRAMP = 0;
        wc = fc / fs;
        [NUM, DEN] = butter(5, wc);
    elseif nargin == 4
        INI_PHASE = varargin{1};
        CARRAMP = 0;
        wc = fc / fs;
        [NUM, DEN] = butter(5, wc);
    elseif (nargin == 5 || nargin == 6)
        INI_PHASE = varargin{1};
        CARRAMP = varargin{2};
        wc = fc / fs;
        [NUM, DEN] = butter(5, wc);
    elseif nargin > 6
        INI_PHASE = varargin{1};
        CARRAMP = varargin{2};
        NUM = varargin{3};
        DEN = varargin{4};
    end
end