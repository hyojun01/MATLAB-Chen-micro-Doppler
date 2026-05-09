%% ================================================================
%  Hovering Quadcopter Micro-Doppler Simulation
%  - Based on Chen-style quadrotor UAV geometry
%  - Body + rotor centers: point scattering model
%  - Rotor blades: line-integral scattering model
%  - Constant RCS assumption
%
%  Coordinate definition:
%    Radar coordinates:      (X, Y, Z), radar at Q = [0;0;0]
%    Reference coordinates:  (X',Y',Z'), origin at quadcopter center O
%    Local coordinates:      (x, y, z), same origin and same orientation
%
%  Hovering case:
%    No translation motion
%    No body rotation motion
%    Only rotor blade rotations exist
% ================================================================
clear; close all; clc;

%% 1. Radar Parameters
c = 3e8;                % speed of light
lambda = 0.0517;        % C-band wavelength
fc = c / lambda;        % carrier frequency

fprintf('Carrier frequency fc = %.3f GHz\n', fc/1e9);

%% 2. Radar-to-quadcopter geometry
O_radar = [50; 0; 20];                                      % quadcopter center in radar coordinates

R0 = norm(O_radar);                                         % distance from radr to O
alpha = atan2(O_radar(2), O_radar(1));                      % azimuth angle
beta = atan2(O_radar(3), hypot(O_radar(1), O_radar(2)));    % elevation angle

alpha_deg = rad2deg(alpha);
beta_deg = rad2deg(beta);

fprintf('R0    = %.3f m\n', R0);
fprintf('alpha = %.3f deg\n', alpha_deg);
fprintf('beta  = %.3f deg\n', beta_deg);

% unit LOS vector from radar to quadcopter center
u_LOS = O_radar / R0;

%% 3. Quadcopter geometry: x-configuration
arm = 0.2;      % rotor center distance in x/y directions

rotor_local = [ arm, -arm,  arm, -arm;
                arm, -arm, -arm,  arm;
                0,    0,    0,    0   ];

NR = 4;         % number of rotors

% positive rotation is CCW when viewed from +z axis
spin_sign = [+1,+1,-1,-1];

%% 4. Rotor blade parameters
NB = 2;                     % number of blades per rotor
blade_root = 0.0;           % blade root radius
blade_tip = 0.07;           % blade tip radius
blade_len = blade_tip - blade_root;

blade_width = 0.025;        % blade width, not explicitly used in 1-D line model

% hovering rotor speed
Omega_rps = 100;            % rotor rotation rate [rev/s] for hovering
omega = 2*pi*Omega_rps;     % [rad/s]

% initial rotor phases
% chen notes that random initial rotation angles can change the period
% pattern.
rng(31);
phi0 = 2*pi*rand(1, NR);    % random initial phase per rotor

%% 5. Constant RCS value
% chen mentions drone body RCS around 0.01 m^2 and blade RCS around 0.001
% m^2
sigma_body = 0;          % body point scatterer RCS
sigma_rotorCenter = 0;   % rotor center / motor-hub point scatterer RCS
sigma_blade = 1e-3;         % total RCS scale per blade

%% 6. Slow-time sampling
fd_max_approx = 2 * omega * blade_tip * cos(beta) / lambda;

fprintf('Approx. max blade-tip micro-Doppler = %.1f Hz\n', fd_max_approx);

fs_slow = 2e6;             % slow-time smapling rate
Tobs = 0.1;                 % observation time
t = 0:1/fs_slow:Tobs-1/fs_slow;
Nt = numel(t);

%% 7. Static point-scatterer returns
R_body = R0;
s_body = sigma_body * exp(-1j*4*pi/lambda * R_body) * ones(1, Nt);

s_rotorCenters = zeros(1, Nt);

for j = 1:NR
    hub_radar = O_radar + rotor_local(:, j);
    R_hub = R0 + u_LOS.' * rotor_local(:,j);

    s_rotorCenters = s_rotorCenters + ...
        sigma_rotorCenter * exp(-1j*4*pi/lambda * R_hub) * ones(1,Nt);
end

%% 8. Rotor blade line-integral returns
s_blades = zeros(1, Nt);            % total return from all rotor blades
s_rotorBlades = zeros(NR, Nt);      % blade return grouped by rotor

for j = 1:NR

    % rotor center projection onto radar LOS
    c_proj = u_LOS.' * rotor_local(:,j);

    % static phase term due to quadcopter center and rotor center
    rotor_static_phase = exp(-1j*4*pi/lambda * (R0 + c_proj));

    % total blade return from j-th rotor
    s_rotor_j = zeros(1, Nt);

    for k = 0:NB-1

        % blade angel
        theta_jk = spin_sign(j)*omega*t + phi0(j) + 2*pi*k/NB;

        % blade direction vector
        e_jk = [cos(theta_jk); sin(theta_jk); zeros(1, Nt)];

        % LOS projection
        gamma_jk = u_LOS.' * e_jk;

        % closed-form line-integral result
        sinc_arg = 4/lambda * blade_len/2 * gamma_jk;

        line_integral_jk = ...
            blade_len .* ...
            exp(-1j*4*pi/lambda * blade_len/2 * gamma_jk) .* ...
            sinc(sinc_arg);

        % blade return 
        s_jk = sigma_blade * rotor_static_phase .* line_integral_jk;

        % sum blades belonging to the same rotor
        s_rotor_j = s_rotor_j + s_jk;
    end

    % store total blade return of each rotor
    s_rotorBlades(j, :) = s_rotor_j;

    % sum all rotor blade returns
    s_blades = s_blades + s_rotor_j;
end

%% 9. Total received baseband signal
s_total = s_body + s_rotorCenters + s_blades;

%% 10. Analytical blade-tip Doppler guide curves
fd_tip = zeros(NR*NB, Nt);
idx = 0;

for j = 1:NR
    for k = 1:NB
        idx = idx + 1;

        theta = spin_sign(j)*omega*t + phi0(j) + 2*pi*(k-1)/NB;

        % Tip velocity in local/reference coordinates
        vx = -blade_tip * spin_sign(j) * omega .* sin(theta);
        vy =  blade_tip * spin_sign(j) * omega .* cos(theta);
        vz = zeros(1, Nt);

        v_tip = [vx; vy; vz];

        Rdot_tip = u_LOS.' * v_tip;
        fd_tip(idx, :) = -(2/lambda) * Rdot_tip;
    end
end

%% 11. Plot: time-domain signal
figure('Name','Received Baseband Signal');
plot(t, real(s_total), 'LineWidth', 1.0); grid on;
xlabel('Time [s]');
ylabel('Real\{s_{total}(t)\}');
title('Total Baseband Return: Body + Rotor Centers + Line-Integral Blades');

%% 12. Plot: micro-Doppler spectrogram
winLen = 4096;
win = hamming(winLen, 'periodic');
noverlap = round(0.90 * winLen);
nfft = 32768;

[S, F, Tspec] = spectrogram(s_total, win, noverlap, nfft, fs_slow, 'centered');

SdB = 20*log10(abs(S) ./ max(abs(S(:))) + eps);

figure('Name','Micro-Doppler Spectrogram');
imagesc(Tspec, F, SdB);
axis xy;
ylim([-3000 3000]);
caxis([-60 0]);
grid on;
colorbar;
xlabel('Time [s]');
ylabel('Doppler frequency [Hz]');
title('Hovering Quadcopter Micro-Doppler Signature');

% hold on;
% plot(t, fd_tip, 'w--', 'LineWidth', 0.5);
% hold off;

%% 13. Plot: power spectrum
NfftSpec = 2^nextpow2(Nt);
Sfreq = fftshift(fft(s_total, NfftSpec));
faxis = (-NfftSpec/2:NfftSpec/2-1) * fs_slow/NfftSpec;

figure('Name','Power Spectrum');
plot(faxis, 20*log10(abs(Sfreq)./max(abs(Sfreq)) + eps), 'LineWidth', 1.0);
grid on;
xlim([-4000 4000]);
ylim([-80 5]);
xlabel('Doppler frequency [Hz]');
ylabel('Normalized magnitude [dB]');
title('Power Spectrum of Hovering Quadcopter Return');

%% 14. Plot: 3-D geometry
figure('Name','3-D Geometry');
hold on; grid on; axis equal;

% Radar
plot3(0, 0, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
text(0, 0, 0, '  Q: Radar', 'FontSize', 10);

% Quadcopter center
plot3(O_radar(1), O_radar(2), O_radar(3), 'ro', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'r');
text(O_radar(1), O_radar(2), O_radar(3), '  O: Quadcopter center', 'FontSize', 10);

% LOS
plot3([0 O_radar(1)], [0 O_radar(2)], [0 O_radar(3)], 'k--', 'LineWidth', 1.0);

% Reference/local coordinate axes at O
axisLen = 0.35;
quiver3(O_radar(1), O_radar(2), O_radar(3), axisLen, 0, 0, ...
    'LineWidth', 1.2, 'MaxHeadSize', 1.0);
text(O_radar(1)+axisLen, O_radar(2), O_radar(3), ' x/X''');

quiver3(O_radar(1), O_radar(2), O_radar(3), 0, axisLen, 0, ...
    'LineWidth', 1.2, 'MaxHeadSize', 1.0);
text(O_radar(1), O_radar(2)+axisLen, O_radar(3), ' y/Y''');

quiver3(O_radar(1), O_radar(2), O_radar(3), 0, 0, axisLen, ...
    'LineWidth', 1.2, 'MaxHeadSize', 1.0);
text(O_radar(1), O_radar(2), O_radar(3)+axisLen, ' z/Z''');

% Rotor centers and blades at t = 0
for j = 1:NR
    hub = O_radar + rotor_local(:,j);

    plot3(hub(1), hub(2), hub(3), 'bo', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'b');

    if spin_sign(j) > 0
        dirText = 'CCW';
    else
        dirText = 'CW';
    end

    text(hub(1), hub(2), hub(3), ...
        sprintf('  R%d (%s)', j, dirText), 'FontSize', 9);

    for k = 1:NB
        theta0 = phi0(j) + 2*pi*(k-1)/NB;
        u0 = [cos(theta0); sin(theta0); 0];

        p_root = hub + blade_root * u0;
        p_tip  = hub + blade_tip  * u0;

        plot3([p_root(1) p_tip(1)], ...
              [p_root(2) p_tip(2)], ...
              [p_root(3) p_tip(3)], ...
              'LineWidth', 2.0);
    end
end

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Radar, Reference/Local Coordinates, Rotor Centers, and Blades at t = 0');
view(35, 25);

%% 15. Summary print
fprintf('\n===== Simulation Summary =====\n');
fprintf('Rotor configuration: x-configuration\n');
fprintf('Rotor centers [m]:\n');
disp(rotor_local);
fprintf('Rotor 1,2: CCW / Rotor 3,4: CW\n');
fprintf('NB = %d blades per rotor\n', NB);
fprintf('Blade root = %.3f m, tip = %.3f m\n', blade_root, blade_tip);
fprintf('Hovering rotor rate = %.1f r/s\n', Omega_rps);
fprintf('Approx. fd max = %.1f Hz\n', fd_max_approx);
fprintf('Slow-time fs = %.1f Hz, Tobs = %.3f s\n', fs_slow, Tobs);