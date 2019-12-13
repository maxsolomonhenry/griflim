# griflim
Spectrogram Inversion with the Griffin and Lim Algorithm

This code features two implementations of the Griffin and Lim algorithm for reconstruction of a signal given a magnitude-only STFT.  The original algorithm is implemented in griffinlim.m, while an updated version for realtime applications (Zhu et al 2007) is implemented in rtisi2.m.  What follows is a brief description of the included scripts:

glimdemo.m --> example application of the griffinlim function

glimwin.m --> creates window whose square overlap-add sums to unity.  For use in LSEE-MSTFT with the Griffin-Lim algorithm.

griffinlim.m --> spectrogram inversion using the Griffin-Lim algorithm (Griffin & Lim 1984).  This function expects a 2x's buffered STFT such as those generated by stft.m (included)

LSEEMSTFT.m --> performs a single inverse STFT using the least squares error estimate (LSEE) algorithm defined by Griffin and Lim (1984).

rtisi2.m --> implementation of the real-time iterative spectrogram inversion algorithm as outlined by Zhu et al (2007).

rtisi2demo.m --> example application of the rtisi2 function

stft.m --> accepts a mono signal x, a user-provided window, and a specified hop size (in samples).  Signal frames are zero-padded (by a factor of 2) and shifted (for frame-centered phase response). 


Thank you for your interest.

// Max Henry

// Department of Music Technology

// McGill University

// Montreal, Qc

// December 13, 2019.
