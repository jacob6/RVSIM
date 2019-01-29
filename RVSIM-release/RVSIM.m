function RVSIM_index = RVSIM(ImageReference, ImageDistored)
% ========================================================================
% RVSIM Index, Version 1.0
% Copyright(c) 2017 Guangyi Yang, Deshi Li, Fan Lu, Yue Liao and Wen Yang
% All Rights Reserved.
%
% ------------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
% ------------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the Riesz 
% transform and Visual contrast sensitivity-based feature SIMilarity
% (RVSIM) index between two images.
%
% Please refer to the following paper
%
% Yang G, Li D, Lu F, et al. RVSIM: a feature similarity method for full-reference image quality assessment[J]. EURASIP Journal on Image and Video Processing, 2018, 2018(1): 6.
%
% Usage:
%    RVSIM_index = RVSIM(img1, img2);
%    
%    The function takes 2 grayscale images img1 and img2 as input. Besides
%    it requires that many values of visual Contrast Sensitivity
%    Function(CSF) to accomplish band weighting. In file
%    "RVSIM_csf-8bands-m2.1min3.mat" we have already stored a set of CSF 
%    values for a log-gabor filter bank whose arguments are shown below.
%
%    nscale          = 5;
%    norient         = 1;
%    minWaveLength   = 3;
%    mult            = 2.1;
%    sigmaOnf        = 0.55;
%    dThetaOnSigma   = 1.5;
%    
%    If any of the arguments(except nscale) should be changed, the band
%    distribution of the filter bank might be changed and hence a new set
%    of CSF values would be necessary in matching with the new filter bank.
%
% ------------------------------------------------------------------------

% CSF = @(omega)2.6*(0.0192+0.114*120*omega).*exp(-(0.114*120*omega).^1.1);
% formula of CSF
load('RVSIM_csf-8bands-m2.1min3.mat')

L = 255;
K1 = 1.09;
K2 = 1.162;
K3 = 1.00;
C1 = (K1*L)^2;
C2 = (K2*L)^2;
C3 = (K3*L)^2;
% Coefficients adjusting the contrast response

filterbankargs.nscale = 5;
filterbankargs.norient = 1;
nb = filterbankargs.nscale*filterbankargs.norient;
[Ib_ref, ~] = loggaborconv2(ImageReference, filterbankargs);
[Ib_dis, ~] = loggaborconv2(ImageDistored, filterbankargs);
% Band-limited images taken through log-gabor filter bank

Rx_ref = cell(1,nb);
Ry_ref = cell(1,nb);
Rx_dis = cell(1,nb);
Ry_dis = cell(1,nb);
for s = 1:nb
    [Rx_ref{s}, Ry_ref{s}] = RieszTransform(Ib_ref{s});
    [Rx_dis{s}, Ry_dis{s}] = RieszTransform(Ib_dis{s});
end
% Riesz transform in 2 directions

A_ref = cell(1,nb);
Ai_ref = cell(1,nb); % imaginary part amplitude
A_dis = cell(1,nb);
Ai_dis = cell(1,nb);
for s = 1:nb
    A_ref{s} = sqrt(Rx_ref{s}.^2 + Ry_ref{s}.^2 + Ib_ref{s}.^2);
    Ai_ref{s} = sqrt(Rx_ref{s}.^2 + Ry_ref{s}.^2);
    A_dis{s} = sqrt(Rx_dis{s}.^2 + Ry_dis{s}.^2 + Ib_dis{s}.^2);
    Ai_dis{s} = sqrt(Rx_dis{s}.^2 + Ry_dis{s}.^2);
end
% Local amplitude of monogenic signals

[G_ref]=gradientxy(ImageReference);
[G_dis]=gradientxy(ImageDistored);
Sim_G = (2 * G_dis .* G_ref + C2)./(G_dis.^2 + G_ref.^2 + C3);
% Gradient magnitude similarity

RVSIM_index = 0;
MPC = MonogenicPC(ImageReference, filterbankargs); % monogenic phase congruence
for s = 1:nb
    Sim_A = (2*A_ref{s}.*A_dis{s} + C1)./(A_ref{s}.*A_ref{s}+A_dis{s}.*A_dis{s} + C1);
    % local amplitude similarity
    Sim_theta = exp(-abs((Rx_ref{s}.*Ry_dis{s} - Ry_ref{s}.*Rx_dis{s})./(Rx_ref{s}.*Rx_dis{s} + Ry_ref{s}.*Ry_dis{s})));
    % local direction similarity
    Sim_phi = exp(-abs((Ai_ref{s}.*Ib_dis{s} - Ib_ref{s}.*Ai_dis{s})./(Ib_ref{s}.*Ib_dis{s} + Ai_ref{s}.*Ai_dis{s})));
    % local phase similarity
    Sim_MS = Sim_theta.*Sim_phi.*Sim_A;
    RVSIMb_index = sum(sum(Sim_MS.*Sim_G.*MPC));
    RVSIM_index = RVSIM_index + csf(s).*RVSIMb_index;
    % band weighting for similarity index
end
RVSIM_index = RVSIM_index./sum(csf);
return;
end