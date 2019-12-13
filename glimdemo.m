% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Sample implementation of griffinlim.m
%
% Accepts a sound specified by 'filename,' performs an STFT using stft.m,
% strips the STFT of its phase, and reconstructs using griffinlim.m
% Plays and saves resulting sound as output.wav, generates plots to compare
% analysis and original signal, and to track the squared errors reduction
% over synthesis iterations.
% 
%
% Based on the algorithm developed by Griffin and Lim:
% D. Griffin and J. Lim. Signal estimation from modified short-time
%     Fourier transform. IEEE Trans. Acoust. Speech Signal Process.,
%     32(2):236-243, 1984.

clearvars;

filename = 'demo.wav';          % input filename here
numIts = 100;                   % specify number of synthesis iterations

[x, fs] = audioread(filename);
x = x(:,1); % take left channel only

X = stft(x, glimwin(4096,4), 4096/4);   % generate stft
X = abs(X);                             % strip STFT of original phase

[y, D] = griffinlim(X, numIts, 4);

sound(y, fs);                           % play synthesis
audiowrite('output.wav', y, fs);        % write to output.wav

% PLOTS
figure();
plot((0:length(x)-1)/fs, x, (0:length(y)-1)/fs, y);
title('GLIM reconstruction');
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