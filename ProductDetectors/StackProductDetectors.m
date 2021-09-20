function [PD1, PD2, PD3, PD4, PD5, PD6, PD1f, PD2f, PD3f, PD4f, PD5f, PD6f, medianSignal, medianSignalShaped] = StackProductDetectors(w, fc, fs, initialPhase, N, doPlot)
% stacks the result of six AM demodulators
% 1. 4th-order trig product detector
% 2. two-copy multiply and shift
% 3. basic product detector
% 4. high-bandwidth lockin demodulator
% 5. matlab builtin (amdemod)
% 6. andreas' envelope detector

    % use highpass filter to detrend signal
    [b, a] = butter(1, 10 / (fs/2), 'high');
    w = filtfilt(b, a, w);

    % use default lowpass filter (5th order butterworth) 
    wc = fc / fs;
    [b, a] = butter(5, wc);

    % trig identity: sin^2(t) * cos^2(t) = (1 = cos(4*t)) / 8
    PD1 = TrigProductDetector4(w(1:N), fc, fs, initialPhase, 0, b, a);

    % two-copy multiply and shift (product detector #2)
    PD2 = MultiplyAndShiftDetector(w(1:N), fc, fs, initialPhase, 0, b, a);

    % basic product detector (product detector #3)
    PD3 = BasicProductDetector(w(1:N), fc, fs, initialPhase, 0, b, a);

    % high bandwidth lockin demodulator (product detector #4)
    PD4 = HighBWLockInDetector(w(1:N), fc, fs, initialPhase, 0, b, a);

    % matlab builtin amdemod (product detector #5)
    try
        PD5 = amdemod(detrend(w(1:N)), fc, fs, initialPhase, 0, b, a);
    catch
        fprintf('Signal_Blocks License Unavailable\n');
    end

    % andreas envelope detector
    PD6 = andreasEnvelopeDetectorIrigB(detrend(w(1:N)), fc, fs);
    
    dt = 0:1/fs:(length(w) / fs - 1 / fs); dt = dt';

    % apply pulse shaper
    PD1f = PadAndShapeFilter(PD1, fs);
    PD2f = PadAndShapeFilter(PD2, fs);
    PD3f = PadAndShapeFilter(PD3, fs);
    PD4f = PadAndShapeFilter(PD4, fs);
    PD6f = PadAndShapeFilter(PD6, fs);
    
    try
        if (doPlot)
            hp = plot(dt(1:N), [PD1, PD2, PD3, PD4, PD6, PD5]);
        end
        medianSignal = median([PD1, PD2, PD3, PD4, PD6, PD5], 2);
        PD5f = PadAndShapeFilter(PD5, fs);
        % medianSignalShaped = median([PD1f,PD2f,PD3f,PD4f,PD5f,PD6f], 2);
    catch
        if (doPlot)
            hp = plot(dt(1:N), [PD1, PD2, PD3, PD4, PD6]);
        end
        medianSignal = median([PD1, PD2, PD3, PD4, PD6], 2);
        % medianSignalShaped = median([PD1f,PD2f,PD3f,PD4f,PD6f], 2);
    end
    medianSignalShaped = PadAndShapeFilter(medianSignal, fs);
    
    if (doPlot)
        hax = gca;
        hax.XLim = [0, 0.11];
        hax.YLim(1) = 0;
        hax.XAxis.Label.String = 'Time (s)';
        hax.YAxis.Label.String = 'Phase (rad)';
        legendStrings = {'sin^2(t) * cos^2(t) = (1 = cos(4*t)) / 8', 'two-copy multiply-and-shift', 'basic product detector', 'high bandwidth lock-in demodulator', 'Andreas-Envelope'};
        if exist('PD5','var')
            legendStrings{end+1} = 'amdemod (MATLAB builtin)';
        end
        legend(hax, legendStrings);
    end

end