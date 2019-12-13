function X = stft(x, window, hop)
% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Generates spectrogram of given input signal 'x', using window length 'window' defined
% by user, and hop size 'hop'.  Each stft frame is generated from a
% windowed, 2x's zero-padded, and shifted signal (for correct phase
% estimation).
%
% This code is heavily indebted to Malcolm Slaney's
% ComplexSpectrogram.m


% Check for proper orientation of window vector
[r, c] = size(window);
if r > c
    window = window';
end

% Likewise for the signal x
[r, c] = size(x);
if r > c
    x = x';
end

windowSize = length(window);
numFrames = floor((length(x) - windowSize)/hop) + 1;

% Calculations for for centered, zero-padded fft frames.
fftLength = 2*windowSize;
fftS = floor(windowSize/2);             % index for start of windowed signal
fftE = fftS + windowSize - 1;           % index for end of windowed signal
fftBuffer = zeros(1, fftLength);

% Initialize output
X = zeros(fftLength, numFrames);

for curFrame = (1:numFrames)
    waveS = 1 + (curFrame - 1)*hop;     % find start/end of signal to be windowed
    waveE = waveS + windowSize - 1;
    
    fftBuffer = 0*fftBuffer;            % clear buffer
    fftBuffer(fftS:fftE) = ...          % place windowed signal in buffer
        x(waveS:waveE) .* window;
    fftBuffer = fftshift(fftBuffer);    % rotate to center

    % write fft frame to output
    X(:, curFrame) = fft(fftBuffer).';
end