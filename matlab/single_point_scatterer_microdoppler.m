%% Motion model:
%   1) Target center translation: R_center(t) = R0 + V*t
%   2) Scatterer rotation around targert center:
%       r_ref(t) = R_rot(t) * R_init * r0_local
%   3) Total scatterer position:
%       p(t) = R_center(t) + r_ref(t)
%   4) Baseband return:
%       s(t) = exp(j * 4*pi/lambda * R(t))
%   5) Doppler:
%       fD(t) = (2/lambda) * n(t)^T * (V + omega_ref x r_ref(t))

clear; close all; clc;

%% 1. Radar parameters
c = 3e8;                        % speed of light
fc = 10e9;                      % carrier frequency, Chen example uses 10GHz
lambda = c / fc;                % wavelength

PRF = 2000;                     % slow-time sampling rate / pulse repetition frequency
N = 4096;                       % number of slow-time samples
t = (0:N-1) / PRF;
dt = 1 / PRF;

SNR_dB = inf;

%% 2. Geometry: radar and target center
% Radar is assumed to be at origin Q = [0,0,0].
R0 = [1000; 5000; 5000];                            % initial target center position

% translation velocity of target center in radar coordinates
n0 = R0 / norm(R0);                                 % far-field initial LOS unit vector
V_trans = 5 * n0;                                   % translation velocity

% if you want arbitrary translation, replace with for example:
% V_trans = [0; -5; 0];

%% 3. Point scatterer in target local coordinates
r0_local = [1.0; 0.6; 0.8];                         % scatterer position in target local coordinates
sigma = 1.0;                                        % scatterer reflectivity amplitude

%% 4. Initial attitude: Chen Euler convention
phi = deg2rad(30);                                  % rotation about z-axis
theta = deg2rad(20);                                % rotation about x-axis
psi = deg2rad(20);                                  % rotation about z-axis

R_init = Rz(phi) * Rx(theta) * Rz(psi);

%% 5. Angular Velocity
% define angular velocity in target local coordinates
omega_local = [pi; 2*pi; pi];                       % [rad/s] in target local coordinates

% convert to reference/radar coordinates
Omega = norm(omega_local);                          % scalar angular speed
omega_axis_ref = R_init * (omega_local / Omega);    % angular velocity vector in reference coordinates
omega_ref = Omega * omega_axis_ref;                 % angular velocity vector in reference coordinates

rotation_period = 2*pi / Omega;

fprintf('Angular speed |omega| = %.4f rad/s\n', Omega);
fprintf('Rotation period T = %.4f s\n', rotation_period);
fprintf('Wavelength lambda = %.4f m\n', lambda);

%% 6. Allocate arrays
r_ref      = zeros(3, N);      % scatterer position relative to target center [m]
p_center   = zeros(3, N);      % target center position [m]
p_scatter  = zeros(3, N);      % absolute scatterer position [m]
R_range    = zeros(1, N);      % radar-to-scatterer range [m]

v_rot      = zeros(3, N);      % rotational velocity omega x r [m/s]
v_total    = zeros(3, N);      % total scatterer velocity [m/s]

fD_exact   = zeros(1, N);      % exact Doppler using instantaneous LOS [Hz]
fD_far     = zeros(1, N);      % far-field Doppler using fixed LOS n0 [Hz]
fD_trans   = zeros(1, N);      % translation Doppler [Hz]
fD_micro   = zeros(1, N);      % rotation-induced micro-Doppler [Hz]

%% 7. Generate scatterer motion and Doppler
for k = 1:N
    tk = t(k);

    % Rodrigues rotation about omega_axis_ref by angle Omega*t
    R_rot = rodrigues(omega_axis_ref, Omega * tk);

    % Scatterer location relative to target center in reference coordinates
    r_ref(:,k) = R_rot * R_init * r0_local;

    % Target center translation
    p_center(:,k) = R0 + V_trans * tk;

    % Absolute scatterer position
    p_scatter(:,k) = p_center(:,k) + r_ref(:,k);

    % Range from radar at origin
    R_range(k) = norm(p_scatter(:,k));

    % Velocity components
    v_rot(:,k) = cross(omega_ref, r_ref(:,k));
    v_total(:,k) = V_trans + v_rot(:,k);

    % Exact LOS changes slightly with time
    n_exact = p_scatter(:,k) / R_range(k);

    % Exact Doppler
    fD_exact(k) = (2/lambda) * dot(v_total(:,k), n_exact);

    % Far-field Chen-style decomposition
    fD_trans(k) = (2/lambda) * dot(V_trans, n0);
    fD_micro(k) = (2/lambda) * dot(v_rot(:,k), n0);
    fD_far(k) = fD_trans(k) + fD_micro(k);
end

%% 8. Generate complex baseband return
phase_total = 4*pi/lambda * (R_range - R_range(1));         % remove phase offset from initial range
s_clean = sigma * exp(1j * phase_total);

% Optional additive complex Gaussian noise
if isfinite(SNR_dB)
    sigPow = mean(abs(s_clean).^2);
    noisePow = sigPow / (10^(SNR_dB/10));
    noise = sqrt(noisePow/2) * (randn(size(s_clean)) + 1j*randn(size(s_clean)));
    s = s_clean + noise;
else
    s = s_clean;
end

% Bulk-translation-compensated signal
R_center = vecnorm(p_center, 2, 1);
phase_center = 4*pi/lambda * (R_center - R_center(1));
s_micro_only = s_clean .* exp(-1j * phase_center);

%% 9. Numerical instantaneous frequency check
phase_unwrapped = unwrap(angle(s_clean));
f_inst_numeric = [diff(phase_unwrapped)/(2*pi*dt), NaN];

%% 10. Plot Doppler curves
figure('Name','Doppler Components');
plot(t, fD_exact, 'LineWidth', 1.2); hold on;
plot(t, fD_far, '--', 'LineWidth', 1.2);
plot(t, fD_trans, ':', 'LineWidth', 1.2);
plot(t, fD_micro, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Single Point Scatterer Doppler Components');
legend('Exact total Doppler', ...
       'Far-field total Doppler', ...
       'Translation Doppler', ...
       'Rotation-induced micro-Doppler', ...
       'Location', 'best');

figure('Name','Numerical Instantaneous Frequency Validation');
plot(t, fD_exact, 'LineWidth', 1.2); hold on;
plot(t, f_inst_numeric, '--', 'LineWidth', 1.0);
grid on;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Doppler from Formula vs. Phase Derivative');
legend('Formula-based exact Doppler', ...
       'Numerical derivative of signal phase', ...
       'Location', 'best');

%% 11. Plot 3D trajectory
figure('Name','3D Scatterer Trajectory');
plot3(p_scatter(1,:), p_scatter(2,:), p_scatter(3,:), 'LineWidth', 1.2); hold on;
plot3(p_center(1,:), p_center(2,:), p_center(3,:), '--', 'LineWidth', 1.2);
grid on; axis equal;
xlabel('U (m)');
ylabel('V (m)');
zlabel('W (m)');
title('3D Trajectory: Target Center and One Rotating Scatterer');
legend('Point scatterer trajectory', 'Target center trajectory', 'Location', 'best');

%% 12. Time-frequency analysis: total Doppler
winLen = 256;
overlap = 240;
nfft = 2048;

[S,F,Tstft] = spectrogram(s, hamming(winLen), overlap, nfft, PRF, 'centered');

figure('Name','STFT: Total Doppler');
imagesc(Tstft, F, 20*log10(abs(S) + eps));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('STFT of Returned Signal: Translation + Micro-Doppler');
colorbar;
hold on;
plot(t, fD_far, 'w', 'LineWidth', 1.5);
ylim([-PRF/2, PRF/2]);

%% Local function
function S = skew3(a)
    % Skew-symmetric matrix such that skew3(a)*b = cross(a,b)
    S = [  0    -a(3)   a(2);
          a(3)   0     -a(1);
         -a(2)  a(1)    0   ];
end

function R = rodrigues(axis_unit, angle_rad)
    % Rodrigues rotation formula:
    % R = I + K sin(theta) + K^2 (1 - cos(theta))
    axis_unit = axis_unit / norm(axis_unit);
    K = skew3(axis_unit);
    R = eye(3) + K*sin(angle_rad) + K*K*(1 - cos(angle_rad));
end

function R = Rz(a)
    R = [cos(a), -sin(a), 0;
         sin(a),  cos(a), 0;
         0,       0,      1];
end

function R = Rx(a)
    R = [1, 0,       0;
         0, cos(a), -sin(a);
         0, sin(a),  cos(a)];
end