function spatial_spec = MUSIC(XCov_FMM, f_range,...
    theta_range, steerVec_base, source_num)
%% -------------------------------------------------------------------
% multiple source classification(MUSIC)
% 
% Usage:
%   spatial_spectrum = MUSIC(XCov_FMM, f_range, theta_range, steerVec)
%
% Inputs:
%   XCov_FMM: spatial corvariance matrix 
%   f_range: slected frequency bins to implement doa estimation
%   theta_range: scan range 
%   steerVec_base: steering vector of the first frequency bin, shape of
%   M*theta_range
%   source_num: number of sources
%
% Output:
%   spatial_spectrum: pedudo spatial spectrum estimated with MUSIC method
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
M = size(XCov_FMM, 3);
%--------------------------------------------------------------------------
%% steer response power
for theta_index = 1: theta_length
    for f_index = 1: f_length
        f = f_range(f_index);
        steerVec = steerVec_base(:, theta_index).^f;
        XCov_MM = squeeze(XCov_FMM(f, :, :));
        [E, D] = eig(XCov_MM);
        [~, index] = sort(diag(D), 'ascend');
        E_sorted = E(:, index);
        Un = E_sorted(:, 1:M-source_num);
        power(theta_index, f_index) =1/ abs(steerVec' * (Un *Un')...
            * steerVec);
    end
end
spatial_spec = mean(power, 2);
spatial_spec = 10 * log10(spatial_spec / max(spatial_spec));
