function y = rtisi2(X, numIts)
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
outputFrames = zeros(fftLength, numFrames);         % output frames, accessed throughout synthesis
fftbuf = zeros(1, fftLength);                       % fft buffer

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
    
    % add tail of prev frame into current frame, and rewindow
    if curFrame > 1
        prevFrame = outputFrames((fftS:fftE), curFrame - 1)';               % window prev frame
        prevFrameContribution = [prevFrame(hop+1:end) zeros(1, hop)];       % extract the tail that overlaps current frame
        withPrevFrameContr = window .* ...                                  % combine with current frame and window
            (prevFrameContribution + outputFrames((fftS:fftE), curFrame)');
        
        outputFrames((fftS:fftE), curFrame) = withPrevFrameContr;           % store overlapped frame back in outputFrames
    end
    
    for i = (1:numIts)
        fftbuf = fftbuf*0;
        
        if curFrame == 1
            % generate random phase for first frame
            newPhasefft = exp(j*2*pi*rand(fftLength, 1))';
        else
            % grab section of the lookahead buffer to recalculate phase
            newphasesource = outputFrames((fftS:fftE),curFrame)';
            fftbuf(fftS:fftE) = newphasesource;
            newPhasefft = fft(fftshift(fftbuf));
        end
        
        % "magnitude-constrained" phase update
        newPhaseMagfft = max(newPhasefft, epsilon) ...
            ./ max(abs(newPhasefft), epsilon) .* goalMag(:, curFrame)';
        
        % place updated frame in outputFrames
        updatedPreWindow = fftshift(real(ifft(newPhaseMagfft)));
        outputFrames((fftS:fftE),curFrame) = updatedPreWindow(fftS:fftE);
        
        % grab signal from center and add to the sum buffer
        newframe = outputFrames((fftS:fftE),curFrame)';
        
        % add tail of prev frame into curFrame and rewindow
        if curFrame > 1
            prevFrame = outputFrames((fftS:fftE), curFrame - 1)';
            prevFrameContribution = [prevFrame(hop+1:end) zeros(1, hop)];
            withPrevFrameContr = window .* ...
                (prevFrameContribution + outputFrames((fftS:fftE), curFrame)');
            
            outputFrames((fftS:fftE), curFrame) = withPrevFrameContr;
        end
    end
    y(wavS:wavE) = y(wavS:wavE) + outputFrames((fftS:fftE),curFrame)';
end
end
