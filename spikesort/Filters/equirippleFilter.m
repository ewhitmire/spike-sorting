function Hd = equirippleFilter
%EQUIRIPPLEFILTER Returns a discrete-time filter object.

%
% MATLAB Code
% Generated by MATLAB(R) 7.14 and the DSP System Toolbox 8.2.
%
% Generated on: 26-Mar-2014 17:02:57
%

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 24000;  % Sampling Frequency

N      = 18;    % Order
Fstop1 = 250;   % First Stopband Frequency
Fpass1 = 300;   % First Passband Frequency
Fpass2 = 2600;  % Second Passband Frequency
Fstop2 = 2800;  % Second Stopband Frequency
Wstop1 = 1;     % First Stopband Weight
Wpass  = 5;     % Passband Weight
Wstop2 = 2;     % Second Stopband Weight
dens   = 30;    % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fstop1 Fpass1 Fpass2 Fstop2 Fs/2]/(Fs/2), [0 0 1 1 0 ...
           0], [Wstop1 Wpass Wstop2], {dens});
Hd = dfilt.dffir(b);

% [EOF]
