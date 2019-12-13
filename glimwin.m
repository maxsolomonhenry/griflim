function window = glimwin(L, OL)
% Max Henry
% Final Project for MUMT 605
% Digital Sound Synthesis and Audio Processing
% Prof. Philippe Depalle
%
% Defines window of length L --as specified in Griffin and Lim (1984)-- 
% that has the property of summing to 1 when squared and overlap-added
% at a rate of OL.
%
% Based on the algorithm developed by Griffin and Lim:
% D. Griffin and J. Lim. Signal estimation from modified short-time
%     Fourier transform. IEEE Trans. Acoust. Speech Signal Process.,
%     32(2):236-243, 1984.

S = floor(L/OL);

wr = ones(1, L) .* (S/L)^(1/2);   % rectangular base

a = 0.54;
b = -0.46;
phi = pi/L;

window = 2*wr/sqrt(4*a^2 + 2*b^2) .* ...            % window forumula
    (a + b * cos(2*pi.*(0:L-1)./L + phi));
end