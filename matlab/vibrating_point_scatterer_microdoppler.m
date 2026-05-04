%% Micro-Doppler Simulation of One vibrating point scatterer in 3D

clear; close all; clc;

%% 1. Radar / slow-time setup
c = 3e8;                % speed of light
fc = 10e9;              % carrier frequency 
lambda = c / fc;        % wavelength

PRF = 2000;             % pulse repetition frequency
Np = 2048;              % number of pulses

Ts = 1 / PRF;           % slow-time sampling interval
t = (0:Np-1) * Ts;      % slow-time vector

%% 2. 3D geometry
% Radar coordinate system: (U, V, M)
% Radar is at origin.
% Target vibration center O is stationary at R0.
R0 = [1000; 5000; 5000];        % target vibration center position

% LOS unit vector from radar to target center
n_LOS = R0 / norm(R0);

r_static = [0; 0; 0];           % static position relative to vibration center

%% 3. Vibration motion setup
% sinusoidal vibration:
%   r_vib(t) = Dv * sin(2*pi*fv*t) * nV
Dv = 0.01;                      % vibration amplitude
fv = 2;                         % vibration frequency

az_vib = deg2rad(20);           % vibration direction azimuth
el_vib = deg2rad(10);           % vibration direction elevation

% Unit vector of vibration direction in 3D reference coordinates
nV = [cos(el_vib)*cos(az_vib); cos(el_vib)*sin(az_vib); sin(el_vib)];

% projection of vibration direction onto radar LOS
proj_LOS = dot(nV, n_LOS);

fprintf('LOS projection nV^T n_LOS = %.6f\n', proj_LOS);

%% 4. Point scatterer position and radial displacement
% Full 3D point scatterer position relative to radar:
%   R_P(t) = R0 + r_static + Dv*sin(2*pi*fv*t)nV
% For micro-Doppler phase, only the scalar range variation matters.
% Under far-field approximation:
%   delta_R(t) ≈ n_LOS^T [Dv*sin(2*pi*fv*t)*nV]
%              = Dv*sin(2*pi*fv*t)*(nV^T n_LOS)
vib_disp = Dv * sin(2*pi*fv*t);         % scalar vibration displacement
delta_R = vib_disp * proj_LOS;          % radial range variation

% full 3D trajectory, useful for debugging/visualization
rP = R0 + r_static + nV * vib_disp;

%% 5. Baseband received signal
% For monostatic radar, round-trip phase is:
%   phi(t) = 4*pi/lambda * R(t)
%
% Static range phase is removed because it does not affect micro-Doppler.
% Therefore:
%   s(t) = exp(j * 4*pi/lambda * delta_R(t))
%
% This is a phase-modulated signal:
%   s(t) = exp(j * beta * sin(2*pi*fv*t))
%
% where:
%   beta = 4*pi/lambda * Dv * (nV^T n_LOS)
sigma = 1.0;                            % point scatterer reflectivity
beta = (4*pi/lambda) * Dv * proj_LOS;   % phase modulation index

s_rx = sigma * exp(1j * beta * sin(2*pi*fv*t));

fprintf('Phase modulation index beta = %.6f rad\n', beta);

%% 6. Theoretical micro-Doppler frequency
% Instantaneous Doppler frequency 
%   f_md(t) = (1/2*pi) * d(phi(t))/dt
%
% Since:
%   phi(t) = 4*pi/lambda * Dv * sin(2*pi*fv*t) * (nV^T n_LOS)
%
% then:
%   f_md(t) = 4*pi*Dv*fv/lambda * cos(2*pi*fv*t) * (nV^T n_LOS)

fmd_theory = (4*pi*Dv*fv/lambda) * proj_LOS .* cos(2*pi*fv*t);

fprintf('Theoretical max |micro-Doppler| = %.3f Hz\n', max(abs(fmd_theory)));

%% 7. Numerical instantaneous frequency check
phi_unwrapped = unwrap(angle(s_rx));
f_inst_num = [diff(phi_unwrapped)/(2*pi*Ts), NaN];

figure;
plot(t, fmd_theory, 'LineWidth', 1.5); hold on;
plot(t, f_inst_num, '--', 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Micro-Doppler Frequency [Hz]');
title('Theoretical vs. Numerical Instantaneous Micro-Doppler');
legend('Theoretical f_{md}(t)', 'Numerical instantaneous frequency');

%% 8. Spectrum: Bessel sideband structure
Nfft = 2^nextpow2(8*Np);
S = fftshift(fft(s_rx, Nfft));
freq = (-Nfft/2:Nfft/2-1) * PRF / Nfft;

figure;
plot(freq, 20*log10(abs(S)/max(abs(S)) + eps), 'LineWidth', 1.2);
grid on;
xlabel('Baseband Frequency [Hz]');
ylabel('Normalized Magnitude [dB]');
title('FFT Spectrum of Vibrating Point Scatterer Return');
xlim([-30 30]);
ylim([-80 5]);

% Theoretical Bessel line amplitudes
kMax = 10;
k = -kMax:kMax;
lineFreq = k * fv;
lineAmp = abs(besselj(k, beta));

figure;
stem(lineFreq, lineAmp, 'filled', 'LineWidth', 1.2);
grid on;
xlabel('Baseband Frequency [Hz]');
ylabel('|J_k(\beta)|');
title('Theoretical Bessel Sideband Amplitudes');
xlim([-kMax*fv-2, kMax*fv+2]);

%% 9. STFT micro-Doppler signature
winLen = 512;
overlap = round(0.90 * winLen);
nfft_stft = 4096;

window = hamming(winLen, 'periodic');

[Sstft, Fstft, Tstft] = spectrogram(s_rx, window, overlap, nfft_stft, PRF, 'centered');

figure;
imagesc(Tstft, Fstft, 20*log10(abs(Sstft) + eps));
axis xy;
xlabel('Time [s]');
ylabel('Frequency [Hz]');
title('STFT Micro-Doppler Signature of Vibrating Point Scatterer');
ylim([-100 100]);
colorbar;
hold on;
plot(t, fmd_theory, 'w', 'LineWidth', 1.5);
legend('Theoretical f_{md}(t)', 'TextColor', 'w');

%% 10. 3D trajectory visualization
figure;
plot3(rP(1,:), rP(2,:), rP(3,:), 'LineWidth', 1.5); hold on;
plot3(R0(1), R0(2), R0(3), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
grid on; axis equal;
xlabel('U [m]');
ylabel('V [m]');
zlabel('W [m]');
title('3D Trajectory of Vibrating Point Scatterer');
legend('Point scatterer trajectory', 'Vibration center');