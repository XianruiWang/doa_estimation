function spatial_spec = cbf(XCov_FMM, f_range,...
    theta_range, steerVec_base)
%% -------------------------------------------------------------------
% steering response power(SRP) with conventional beamforming(CBF)
% 
% Usage:
%   spatial_spectrum = cbf(XCov_FMM, f_range, theta_range, steerVec)
%
% Inputs:
%   XCov_FMM: spatial corvariance matrix 
%   f_range: slected frequency bins to implement doa estimation
%   theta_range: scan range 
%   steerVec_base: steering vector of the first frequency bin, shape of
%   M*theta_range
%
% Output:
%   spatial_spectrum: spatial spectrum estimated with CBF method
%
% Author:
%   Xianrui Wang, Center of Intelligent Acoustics and Immersive
%   Communications(CIAIC)
%
% Contact:
%   wangxianrui@mail.nwpu.edu.cn
%--------------------------------------------------------------------------
theta_length = length(theta_range);
f_length= length(f_range);
%--------------------------------------------------------------------------
%% steer response power
for theta_index = 1: theta_length
    for f_index = 1: f_length
        f = f_range(f_index);
        steerVec = steerVec_base(:, theta_index).^f;
        XCov_MM = squeeze(XCov_FMM(f, :, :));
        power(theta_index, f_index) = abs(steerVec' * XCov_MM * steerVec);
    end
end
spatial_spec = mean(power, 2);
spatial_spec = 10 * log10(spatial_spec / max(spatial_spec));
