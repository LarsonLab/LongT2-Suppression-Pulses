function [rf, rf_g] = design_dualband_pulse(N, pw, bw, foff2, flips, d);
% DESIGN_DUALBAND_PULSE Creates a maximum-phase dual-band RF suppression 
%        pulse using the complex Parks-McClellan filter design algorithm
%    [rf, rf_g] = design_dualband_pulse(N, pw, bw, foff, [flips, d])
%
%    Uses the complex Parks-McClellan filter design algorithm (see CFIRPM) 
%    to design a dual-band RF pulse.  The input bandwidths are the exact
%    suppression bandwidths,  not the full-width, half-maximum bandwidths.
%    The resulting pulse is maximum phase.
%
%    A filter order of N > 70 is recommended.
% Inputs:
%    N - number of points in pulse / filter order
%    pw - pulsewidth in seconds
%    bw - Bandwidths of suppression bands
%    foff2 - frequency offsets of suppression band, in Hz. -220 for 1.5T water/fat sat
%    flips (optional) - flips of suppression bands (default is [pi/2 pi/2])
%    d (optional) - ripples in suppression and no sat bands (default is [.01 .01])
% Outputs:
%    rf - normalized rf pulse (sums to flip angle)
%    rf_g - rf pulse in Gauss
%
% Example:
%   [rf rf_g] = 
%
% NOTE: fmp.m, b2rf.m, b2a.m, ab2rf.m and rfscaleg.m must be installed to 
% run this function.
%
% Written by Peder Larson, 12/13/2005
% (c) Board of Trustees, Leland Stanford Junior University

if (nargin < 5) || isempty(flips)% no flips supplied
    flips = [pi/2 pi/2];
end

if (nargin < 6) || isempty(d) % no ripple supplied
    d = [.01 .01];
end

% Adjust ripples for suppression pulse (will not be exact when flips are
% not pi/2...)
d1 = d(1)/2;
d2 = sqrt(d(2));

% normalized frequencies
bw1n = bw(1)/2*pw;
bw2n = bw(2)/2*pw;
foff2n = foff2*pw;

% for Max/Min phase pulse
d1m = 2*d1; d2m = d2^2/2;
di = 0.5*dinf(d1m,d2m);

% transition width
bw_trans = di/pw;
% normalized transition width
bw_t = bw_trans*pw;

max_flip = max(flips);

if (foff2n < 0)
  f = [-N/2 foff2n-bw2n-bw_t foff2n-bw2n foff2n+bw2n foff2n+bw2n+bw_t ...
       -bw1n-bw_t -bw1n  bw1n  bw1n+bw_t N/2]/(N/2);
a = [0 0 sin(flips(2)/2) sin(flips(2)/2) 0 0 sin(flips(1)/2), sin(flips(1)/2) 0 0] / sin(max_flip/2);
else
  f = [-N/2 -bw1n-bw_t -bw1n  bw1n  bw1n+bw_t foff2n-bw2n-bw_t foff2n-bw2n ...
       foff2n+bw2n foff2n+bw2n+bw_t N/2]/(N/2);
a = [0 0 sin(flips(1)/2) sin(flips(1)/2) 0 0 sin(flips(2)/2), sin(flips(2)/2) 0 0] / sin(max_flip/2);
end

if any(diff(f) <= 0),
   error(['Pulse with given parameters is not feasible.\n' ...
          'Try decreasing the bandwidths, increasing the pulse duration\n' ...
          'or increasing the ripple values'],1)
end


% Creates linear phase filter
Wpass = d2m/d1m;
w = [1  Wpass 1 Wpass 1];
h_lin =  cfirpm((N-1)*2,f, a, w);
% Converts to maximum phase filter - this compromises the frequency
% response for complex filters!
h = fmp(h_lin(end:-1:1)); % h is of length N

rf = b2rf(h*sin(max_flip/2));
% scale to max_flip
rf = flips(1)*rf/sum(rf);
rf_g = rfscaleg(rf, pw*1e3);
