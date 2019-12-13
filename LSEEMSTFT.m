function x = LSEEMSTFT(X, window, hop, iter, initrand, normalize, debugplots)
% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Reconstructs a signal 'x' from a given STFT representation 'X' and hop
% size 'hop'.
%
% As part of the LSEE algorithm, signal is re-windowed in synthesis then
% normalized by a sum of squared windows over time.  However, the use of a
% self-normalizing window is assumed (see glimwin.m).  If the user would
% like to provide their own window, set 'normalize' to a non-zero
% value.
%
% Based on the algorithm developed by Griffin and Lim:
% D. Griffin and J. Lim. Signal estimation from modified short-time
%     Fourier transform. IEEE Trans. Acoust. Speech Signal Process.,
%     32(2):236-243, 1984.

% This code is heavily indebted to Malcolm Slaney's
% InvertSpectrogram.m

if nargin < 4
    iter = 0;
end
if nargin < 5
    initrand = 1;
end
if nargin < 6
    normalize = 0;
end
if nargin < 7
    debugplots = 0;
end

[r c] = size(window);
if r > c
    window = window';
end

[fftLength, numFrames] = size(X);
windowSize = floor(fftLength/2);

% pointers for moving through zero-padded fft
fftS = floor(windowSize/2);
fftE = fftS + windowSize - 1;

% initialize output signal
x = zeros(1, numFrames * hop + windowSize - 1);

if normalize ~= 0
    % to sum windows (for use in LSEE-MSTFT denominator)
    windowsum = x;
end

for curFrame = (1:numFrames)
    waveS = (curFrame - 1)*hop + 1;
    waveE = waveS + windowSize - 1;
    
    % if selected, initialize first iteration with random phase
    if iter == 1 && initrand == 1        
        X(:, curFrame) = X(:, curFrame)...
            .* exp(j*2*pi*rand(fftLength, 1));
    end
    
    newframe = ifft(X(:, curFrame))';
    shiftframe = real(fftshift(newframe));

    % build up output signal by accumulating windowed frames
    x(waveS:waveE) = x(waveS:waveE) + window.* shiftframe(fftS:fftE);
    
    % plots for debugging/presentation
    if debugplots == 1
        subplot (3, 1, 2); plot(shiftframe); title('ifft shifted'); ...
            xlabel('time index'); ylabel('amplitude'); xlim([1, fftLength]);
        subplot(3, 1, 3); plot(window.* shiftframe(fftS:fftE)); title('windowed'); ...
            xlabel('time index'); ylabel('amplitude'); xlim([1, windowSize]); pause();
        subplot(3, 1, 1); plot(real(newframe)); title('ifft'); ...
            xlabel('time index'); ylabel('amplitude'); xlim([1, fftLength]);
    end
    
    if normalize ~= 0
        % accumulate squared windows for denominator
        windowsum(waveS:waveE) = windowsum(waveS:waveE) + window.^2;
    end
end

if normalize ~= 0
    x = x./windowsum; % normalize by squared window as per LSEE-MSTFT)
end
end