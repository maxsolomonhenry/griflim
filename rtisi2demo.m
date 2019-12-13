% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Sample implementation of rtisi2.m
%
% Accepts a sound specified by 'filename,' performs an STFT using stft.m,
% strips the STFT of its phase, and reconstructs using rtisi2.m
% Generates plots to compare estimate and original signal.  The second half
% of this script may be used to demonstrate how estimate error reduces with
% iterations.  Be careful, as it takes a long time to run.
% 
%
% Based on the algorithm developed by Griffin and Lim:
% D. Griffin and J. Lim. Signal estimation from modified short-time
%     Fourier transform. IEEE Trans. Acoust. Speech Signal Process.,
%     32(2):236-243, 1984.
%
% X. Zhu, G. T. Beauregard, and L. L. Wyse, Real-Time Signal Estimation
%     From Modified Short-Time Fourier Transform Magnitude Spectra.
%     IEEE Trans. Audio Speech Lang. Process.,
%     15(5):1645?1653, 2007.

clearvars;

filename = 'demo.wav';          % input filename here
numIts = 5;                     % specify number of synthesis iterations

% generate a practice STFT
window = glimwin(4096,4);

[x, fs] = audioread(filename);
x = x(:,1); % take left channel only

X = stft(x, window, 4096/4);            % generate stft
X = abs(X);                             % strip STFT of original phase

y = rtisi2(X, numIts);

% measure estimation error D
Xhat = stft(y, window, 4096/4);         % stft of signal estimate
D = sum(sum(abs(X - abs(Xhat))));

% PLOTS
plot((0:length(x)-1)/fs, x, (0:length(y)-1)/fs, y);
title(['RTISI reconstruction, D = ' num2str(D)]);
xlabel('time (s)');
ylabel('amplitude');
legend({'original', 'synthesis'});

sound(y, fs);

%% Test performance over multiple iterations

% Measure estimation error over multiple iterations and generates plots.
% Takes a long time to run!

maxIts = 100;           % specify maximum number of iterations
D = zeros(1, maxIts);   % array to store distance measures

for numIts = 1:maxIts
    y = rtisi2(X, numIts);
    Xhat = stft(y, window, 4096/4);
    D(numIts) = sum(sum(abs(X - abs(Xhat))));
end

plot((0:length(x)-1)/fs, x, (0:length(y)-1)/fs, y);
title('RTISI reconstruction');
xlabel('time (s)');
ylabel('amplitude');
legend({'original', 'synthesis'});

figure();
subplot(2, 1, 1);
plot((0:length(y)-1)/fs, y);
title('Reconstructed signal y[n]');
xlabel('time (s)');
ylabel('amplitude');
subplot(2, 1, 2);
plot(D);
title('Distance measure D over number of iterations');
xlabel('number of iterations');
ylabel('D');