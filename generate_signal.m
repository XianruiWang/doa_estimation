function [multichannel_signal, fs] = generate_signal(room_size, t60,...
    source_pos, array_pos, signal_path)
%% -------------------------------------------------------------------
% Generate multichannel  signal with given room size, reverbertation 
% time(t60), source position and distribution of microphones using toolbox
% RIR-Generator.
%
% Usage:
% multichannel_signal = generate_signal(room_size, t60,...
%                                       source_pos, array_pos, signal_path)
%
% Inputs:
%   room_size: 1x3 array, shape of room
%   t60: reverberation time T60
%   source_pos: 1x3 array, position of source signal
%   array: Mx3 array, position of microphones
%   signal_path: path of input audio signal
%
% Output:
%   multichannel_signal: multichannel signal
%   fs: sampling rate, equal to that of input file
%
% Author:
%   Xianrui Wang, Center of Intelligent Acoustics and Immersive
%   Communications(CIAIC)
%
% Contact:
%   wangxianrui@mail.nwpu.edu.cn
%--------------------------------------------------------------------------
addpath('./RIR-Generator-master');
%--------------------------------------------------------------------------
%% basic configuration
c = 340;                     % sound velocity
nsample = 4096;              % length of impulse response
%--------------------------------------------------------------------------
%% load signal
[signal, fs] = audioread(signal_path);
%% generate impulse response and generate signal
h = rir_generator(c, fs, array_pos, source_pos, room_size, t60, nsample);
M = size(array_pos, 1);
for m = 1:M
    tmp = filter(h(m, :), 1, signal);
    if m == 1
        multichannel_signal = zeros(length(tmp), M);
    end
    multichannel_signal(:, m) = tmp;
end
