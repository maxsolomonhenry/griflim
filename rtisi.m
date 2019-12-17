function y = rtisi(X, numIts)
% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Implementation of the Griffin and Lim algorithm for phase reconstruction
% with Real-Time implementation (RTISI) update by Zhu et. al.
%
% X = input spectrogram (magnitude only STFT) - each frame is assumed to
% be 2-times buffered (zero padded) and centered.  The included script
% stft.m may be used to accomplish this.
%
% numIts = number of synthesis iterations per frame.  10 is usually more
% than sufficient.
%
% D. Griffin and J. Lim. Signal estimation from modified short-time
%     Fourier transform. IEEE Trans. Acoust. Speech Signal Process.,
%     32(2):236-243, 1984.
%
% X. Zhu, G. T. Beauregard, and L. L. Wyse, Real-Time Signal Estimation
%     From Modified Short-Time Fourier Transform Magnitude Spectra.
%     IEEE Trans. Audio Speech Lang. Process.,
%     15(5):1645?1653, 2007.

if nargin < 2
    numIts = 10;
    disp('Default iterations: 10');
end

goalMag = abs(X);                         % force magnitude-only (in case of complex input)
[fftLength, numFrames] = size(goalMag);   % get spectrogram dimensions

% assuming 1/4 overlap (necessary for squared-window-sum normalization assumption)
OL = 4;
windowSize = floor(fftLength/2);
hop = floor(windowSize/OL);

% generate synthesis/analysis window
window = glimwin(windowSize, OL);

% initialize output/calculation buffers
y = zeros(1, numFrames*hop + windowSize - 1);       % output
prevFrame = zeros(1, windowSize);
newFrame = zeros(1, windowSize);
prevFrameBuf = zeros(1, fftLength);
newFrameBuf = zeros(1, fftLength);

% init/define pointers
fftS = floor(windowSize/2);
fftE = fftS + windowSize - 1;
wavS = 0;
wavE = 0;

% a small number (used to avoid division by zero)
epsilon = 1/1e8;

for curFrame = 1:numFrames
    % update output pointers
    wavS = hop * (curFrame - 1) + 1;
    wavE = wavS + windowSize - 1;
    
    % prepare contribution of previous frames
    prevFrame = y(wavS:wavE) .* window;
    prevFrameBuf(fftS:fftE) = prevFrame;
    
    for i = (1:numIts)
        % generate spectrum of what is currently in output signal frame
        prevFrameSpectrum = fft(fftshift(prevFrameBuf));
        
        % replace magnitude of this spectrum with goal magnitude
        updateSpectrum = goalMag(:,curFrame).' .* max(prevFrameSpectrum, epsilon) ./ max(abs(prevFrameSpectrum), epsilon);
        
        %ifft to retrieve new time signal
        newFrameBuf = fftshift(ifft(updateSpectrum));
        newFrame = newFrameBuf(fftS:fftE);
        
        % add to prevFrameBuf in case of iteration
        prevFrameBuf = prevFrameBuf * 0;
        prevFrameBuf(fftS:fftE) = window .* (prevFrame + newFrame);
    end
    
    % overlap add to output signal
    y(wavS:wavE) = y(wavS:wavE) + newFrame;
end
end
