%% Chen-style 3D Rotation-Induced micro-Doppler Simulation
% Case: No translation motion, rigid body rotation only
% Target: 3 point scatterers
%
% Sign convention:
%   This script uses Chen's phase convention:
%       s(t) = exp(+j*4*pi*R(t)/lambda)
%   If your receiver model uses exp(-j*4*pi*R(t)/lambda),
%   the Doppler sign will be flipped.

clear; close all; clc;

%% 1. Radar parameters
c = 3e8;                % speed of light
fc = 10e9;              % carrier frequency
lambda = c / fc;        % wavelength

%% 2. Observation parameters
Tobs = 2.048;           % observation time
N = 8192;               % number of slow-time samples
fs = N / Tobs;          % slow-time sampling rate
t = (0:N-1) / fs;       % time vector

%% 3. Target center / radar LOS
% radar is assumed to be at the origin of radar coordinates.
% target center is fixed: no translation velocity.
target_center = [1000; 5000; 5000];         % [X; Y; Z] in radar coordinates

R0 = norm(target_center);
n_LOS = target_center / R0;                      % unit LOS vector from radar to target center

az = atan2(target_center(2), target_center(1));
el = asin(target_center(3) / R0);

fprintf('LOS azimuth   = %.3f deg\n', rad2deg(az));
fprintf('LOS elevation = %.3f deg\n', rad2deg(el));
fprintf('Initial range = %.3f m\n', R0);

%% 4. Chen initial Euler rotation
% Chen convention:
%   R_init = Rz(phi) * Rx(theta) * Rz(psi)
%
% This is a z-x-z Euler sequence.
phi = deg2rad(30);          % rotation about z-axis
theta = deg2rad(30);        % rotation about x-axis
psi = deg2rad(20);          % rotation about z-axis again

R_init = Rz(phi) * Rx(theta) * Rz(psi);

%% 5. Point scatterer model
% Columns are P0, P1, P2 in local coordinates.
P_local = [...
    0.0, 0.5, -0.3;
    0.0, 0.3, 0.0;
    0.0, 0.4, 0.6];

num_scatterers = size(P_local, 2);

% reflectivity of each point scatterer.
rho = [1.0, 1.0, 1.0];

%% 6. Rigid-body spinning motion
num_rotations = 8;
Omega = 2*pi*num_rotations / Tobs;      % angular speed

% spinning axis in the local coordinates.
spin_axis_local = [0; 0; 1];

% transform local rotation axis into reference coordinates.
spin_axis_ref = R_init * spin_axis_local;
spin_axis_ref = spin_axis_ref / norm(spin_axis_ref);

K = skew3(spin_axis_ref);

fprintf('Omega = %.6f rad/s\n', Omega);
fprintf('Equivalent spin rate = %.6f Hz\n', Omega/(2*pi));
fprintf('Reference rotation axis = [%.4f %.4f %.4f]^T\n', ...
    spin_axis_ref(1), spin_axis_ref(2), spin_axis_ref(3));

%% 7. simulation arrays
pos_ref = zeros(3, num_scatterers, N);      % relative position from target center
vel_ref = zeros(3, num_scatterers, N);      % relative velocity due to rotation
f_mD = zeros(num_scatterers, N);            % theoretical micro-Doppler

s_bb = zeros(1, N);                         % complex baseband return

% To match Chen's far-field formula, use far-field phase:
%   R_k(t) ≈ R0 + n_LOS^T r_k(t)
% Set this false if you want exact spherical range:
%   R_k(t) = || target_center + r_k(t) ||
use_far_field_phase = true;

%% 8. Main rotation loop
for ii = 1:N
    ti = t(ii);

    % rodrigues rotation matrix in reference coordinates
    R_t = eye(3) + K*sin(Omega*ti) + (K*K)*(1 - cos(Omega*ti));

    for kk = 1:num_scatterers
        % position of kth point scatterer relative to target center
        r_k = R_t * R_init * P_local(:, kk);

        % rigid-body rotational velocity
        v_k = Omega * (K * r_k);

        pos_ref(:, kk, ii) = r_k;
        vel_ref(:, kk, ii) = v_k;

        % chen-style roation-induced micro-Doppler
        f_mD(kk, ii) = (2/lambda) * dot(v_k, n_LOS);

        % Range used for phase modulation
        % It is very important for our Project
        if use_far_field_phase
            R_k = R0 + dot(n_LOS, r_k);
        else
            R_k = norm(target_center + r_k);
        end
         % complex baseband return from kth point scatterer
         s_bb(ii) = s_bb(ii) + rho(kk) * exp(1j * 4*pi/lambda * R_k);
    end
end

% Remove strong static/DC component to make micro-Doppler lines easier to see.
% This does not change the theoretical f_mD curves.
s_plot = s_bb;
% s_plot = s_bb - mean(s_bb);

%% 9. Plot theoretical micro-Doppler curves
figure('Name', 'Theoretical micro-Doppler curves');
plot(t, f_mD(1,:), 'LineWidth', 1.2); hold on;
plot(t, f_mD(2,:), 'LineWidth', 1.2);
plot(t, f_mD(3,:), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Micro-Doppler frequency [Hz]');
title('Chen-style rotation-induced micro-Doppler of 3 point scatterers');
legend('P0', 'P1', 'P2', 'Location', 'best');

%% 10. Time-frequency analysis using STFT / spectrogram
win_len  = 256;
noverlap = 240;
nfft     = 2048;

window = kaiser(win_len, 8);

[S, F, Tstft] = spectrogram(s_plot, window, noverlap, nfft, fs, 'centered');

S_dB = 20*log10(abs(S) / max(abs(S(:))) + eps);

figure('Name', 'Micro-Doppler spectrogram');
imagesc(Tstft, F, S_dB);
axis xy;
ylim([-1500 1500]);
clim([-60 0]);
colorbar;
xlabel('Time [s]');
ylabel('Doppler frequency [Hz]');
title('Simulated micro-Doppler signature: rigid-body rotation only');
hold on;

% Overlay theoretical curves
plot(t, f_mD(1,:), 'w--', 'LineWidth', 1.0);
plot(t, f_mD(2,:), 'w--', 'LineWidth', 1.0);
plot(t, f_mD(3,:), 'w--', 'LineWidth', 1.0);

%% 11. 3D scatterer trajectory visualization
figure('Name', '3-D trajectories in reference coordinates');
hold on; grid on; axis equal;

for kk = 1:num_scatterers
    x = squeeze(pos_ref(1, kk, :));
    y = squeeze(pos_ref(2, kk, :));
    z = squeeze(pos_ref(3, kk, :));
    plot3(x, y, z, 'LineWidth', 1.2);
    scatter3(x(1), y(1), z(1), 50, 'filled');
end

xlabel('X relative [m]');
ylabel('Y relative [m]');
zlabel('Z relative [m]');
title('Rotational trajectories of 3 point scatterers');
legend('P0 trajectory', 'P0 start', ...
       'P1 trajectory', 'P1 start', ...
       'P2 trajectory', 'P2 start', ...
       'Location', 'best');

%% Local functions
function R = Rz(a)
    R = [ cos(a), -sin(a), 0;
          sin(a),  cos(a), 0;
               0,       0, 1 ];
end

function R = Rx(a)
    R = [ 1,      0,       0;
          0, cos(a), -sin(a);
          0, sin(a),  cos(a) ];
end

function K = skew3(e)
    % Skew-symmetric matrix such that K*x = e cross x
    e = e(:);
    K = [    0, -e(3),  e(2);
          e(3),     0, -e(1);
         -e(2),  e(1),     0 ];
end