clear;clc;
c = 340;               % sound velocity
L = [7 6 4];           % size of home
%# construct a uniform linear array(ULA)
center_x = 3.5;                     % center of ULA in x axis
center_y = 3;                       % center of ULA in y axis
d = 0.05;                           % interelement spacing
M = 16;                          % number of microphones
left_x = center_x - round(M/2-1)*d; % left end of array
r_x = (left_x : d: left_x+(M-1)*d)';% distribution of microphones on x axis
r_y = repmat(center_y, M, 1);       % location of microphones on y axis
r_z = 2*ones(M, 1);                 % location of microphones on z axis
r = [r_x r_y r_z];                  % location of microphone array
%# location of source signal
doa = 70/180*pi;                    % true direction of arrival
distance = 2;                       % distance between source and center
s_x = center_x + distance*cos(doa);
s_y = center_y + distance*sin(doa);
s_z = 2;
s = [s_x s_y s_z];
%# reverberation time
t60 = 0.4;
%# generate multichannel signal
[multichannel_signal, fs] = generate_signal(L, t60,...
    s, r, 'ref.wav');
%# implement short time fourier transform(STFT)
nFFT = 2^10;
nShift = nFFT/2;
[X_FTM, win] = STFT(multichannel_signal, nFFT, nShift, 'hamming');
%# estimate spatial corvariance matrix(SCM)
[fBin, nFrm, nCh] = size(X_FTM);
XCov_FTMM = X_FTM .* conj(reshape...
    (X_FTM, fBin, nFrm, 1, nCh));    % SCM
XCov_FMM = squeeze(mean(XCov_FTMM, 2));
%# construct steering vector at each theta
f_range = 30: 1: 200;
f_length = length(f_range);
theta_range = 0: 1: 180;
theta_length = length(theta_range);
for theta_index = 1: theta_length
    theta = theta_range(theta_index);
    theta_omega = theta / 180 *pi;
    tau = d * cos(theta_omega) / c;
    f_base = 1 / nFFT * fs * 2* pi;
    steerVec_base(:, theta_index) = exp(1i*f_base*tau.*(0:(M-1))');
end
%# spatial spectrum estimated with CBF method
X_FTM_nor = X_FTM ./ abs(X_FTM);
XCov_FTMM_nor = X_FTM_nor .* conj(reshape...
    (X_FTM_nor, fBin, nFrm, 1, nCh));    % SCM
XCov_FMM_nor = squeeze(mean(XCov_FTMM_nor, 2));
spatial_spec = Capon(XCov_FMM, f_range, theta_range, steerVec_base);
plot(theta_range.', spatial_spec);
[~, doaHat] = max(spatial_spec);









