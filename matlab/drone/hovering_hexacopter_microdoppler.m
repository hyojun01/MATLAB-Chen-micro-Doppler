%% ================================================================
%  Hovering Hexacopter Micro-Doppler Simulation
%  - Body + rotor hubs: point scattering model
%  - Rotor blades: closed-form line-integral scattering model
%  - Far-field approximation for body, hubs, and blades
%
%  Coordinate definition:
%    Radar coordinates:      (X, Y, Z), radar at Q = [0;0;0]
%    Reference coordinates:  (X',Y',Z'), origin at hexacopter center O
%    Local coordinates:      (x, y, z), same origin and same orientation
%
%  Hovering case:
%    No translation motion
%    No body rotation motion
%    Only rotor blade rotations exist
%
%  Important convention:
%    +spin_sign = CCW when viewed from +z axis
%    -spin_sign = CW  when viewed from +z axis
% ================================================================

clear; close all; clc;

%% 1. Radar parameters
c = 3e8;                % speed of light
lambda = 0.0517;        % wavelength, C-band
fc = c / lambda;        % carrier frequency

fprintf('Carrier frequency fc = %.3f GHz\n', fc/1e9);

%% 2. Radar-to-hexacopter geometry
% hexacopter center O expressed in radar coordinates
O_radar = [50; 0; 20];

R0 = norm(O_radar);                                      % distance from radar Q to hexacotper center O
alpha = atan2(O_radar(2), O_radar(1));                   % azimuth angle
beta = atan2(O_radar(3), hypot(O_radar(1), O_radar(2))); % elevation angle

alpha_deg = rad2deg(alpha);
beta_deg  = rad2deg(beta);

fprintf('R0    = %.3f m\n', R0);
fprintf('alpha = %.3f deg\n', alpha_deg);
fprintf('beta  = %.3f deg\n', beta_deg);

% unit LOS vector from radar Q from hexacopter center O
u_LOS = O_radar / R0;

%% 3. Hexacopter geometry: regular hexagon configuration
NR = 6;             % number of rotors

% rotor center distance from hexacopter center O
rotor_radius = sqrt(0.2^2 + 0.2^2);

% hexacopter yaw offset in local coordinates
psi_frame = deg2rad(0);

% rotor angular positions around O in x-y plane
psi_rotor = psi_frame + (0:NR-1) * 2*pi/NR;

% rotor hub positions in local coordinates
rotor_local = rotor_radius * ...
    [cos(psi_rotor);
     sin(psi_rotor);
     zeros(1,NR)];

% alternating spin directions for hexacopter
spin_sign = [+1,-1,+1,-1,+1,-1];

fprintf('\nRotor centers in local/reference coordinates [m]:\n');
disp(rotor_local);

%% 4. Rotor blade parameters
NB = 2;                             % number of blades per rotor

blade_root = 0.0;                   % blade root radius from rotor hub
blade_tip = 0.07;                   % blade tip radius from rotor hub
blade_len = blade_tip - blade_root; % length of blade
blade_mid = 0.5 * (blade_root + blade_tip);

blade_width = 0.025;                % not explicitly used in 1-D line model

% hovering rotor speed
Omega_rps = 100;                    % rotor rotation rate
omega = 2*pi*Omega_rps;             % rotor angluar speed

% initial rotor phases
rng(24);
phi0 = 2*pi*rand(1, NR);

%% 5. Scattering amplitude parameters
sigma_body = 1e-2;                  % body point-scatterer amplitude scale
sigma_rotorHub = 1e-3;              % each rotor hub point-scatterer ampltiude scale
sigma_blade = 2e-3;                 % each blade line-scatterer amplitude scale

removeStaticMeanForTF = true;       % recommended when body/hub returns are nonzero

%% 6. Slow-time sampling
u_xy_norm = norm(u_LOS(1:2));

fd_max_approx = 2 * omega * blade_tip * u_xy_norm / lambda;

fprintf('Approx. max blade-tip micro-Doppler = %.1f Hz\n', fd_max_approx);

fs_slow = 2e6;                      % slow-time sampling rate
Tobs = 0.1;                         % observation time

t = 0:1/fs_slow:Tobs-1/fs_slow;
Nt = numel(t);

fprintf('Slow-time fs = %.1f Hz, Tobs = %.3f s, Nt = %d\n', ...
    fs_slow, Tobs, Nt);

%% 7. static point-scatterer returns: body and rotor hubs
% body
R_body = R0;
s_body = sigma_body * exp(-1j*4*pi/lambda * R_body) * ones(1,Nt);

% rotor hubs
s_rotorHubs = zeros(1,Nt);

for j=1:NR
    R_hub = R0 + u_LOS.' * rotor_local(:,j);

    s_rotorHubs = s_rotorHubs + ...
        sigma_rotorHub * exp(-1j*4*pi/lambda * R_hub) * ones(1,Nt);
end

%% 8. Rotor blade line-integral returns
s_blades = zeros(1,Nt);             % total blade return from all rotors
s_rotorBlades = zeros(NR,Nt);       % blade return grouped by rotor
s_eachBlade = zeros(NR,NB,Nt);      % optional storage for debugging

for j = 1:NR

    % rotor hub projection onto radar LOS
    c_proj = u_LOS.' * rotor_local(:,j);

    % static two-way phase due to center range and rotor hub offset
    rotor_static_phase = exp(-1j*4*pi/lambda * (R0 + c_proj));

    % total blade return from j-th rotor
    s_rotor_j = zeros(1,Nt);

    for k = 0:NB-1

        % blade angle in local coordinates
        theta_jk = spin_sign(j)*omega*t + phi0(j) + 2*pi*k/NB;

        % unit vector along blade direction in lcoal coordinates
        e_jk = [cos(theta_jk);
                sin(theta_jk);
                zeros(1,Nt)];

        % LOS projection of blade direction
        gamma_jk = u_LOS.' * e_jk;

        % closed-form line-integral result
        % closed form:
        %   L * exp(-j*K*r_mid*gamma) *
        %   sinc( K*L*gamma/(2*pi) )
        %
        % MATLAB sinc(x) = sin(pi*x)/(pi*x), so:
        %   K*L/(2*pi) = 2*L/lambda
        sinc_arg = 4/lambda * blade_len/2 .* gamma_jk;

        line_integral_jk = ...
            blade_len .* ...
            exp(-1j*4*pi/lambda * blade_mid .* gamma_jk) .* ...
            sinc(sinc_arg);

        % blade return 
        s_jk = sigma_blade * rotor_static_phase .* line_integral_jk;

        % store individual blade return
        s_eachBlade(j, k+1, :) = reshape(s_jk, 1, 1, []);

        % sum blades belonging to same rotor
        s_rotor_j = s_rotor_j + s_jk;
    end

    % store total blade return from j-th rotor
    s_rotorBlades(j, :) = s_rotor_j;

    % sum all rotor blade returns
    s_blades = s_blades + s_rotor_j;
end

%% 9. Total received baseband signal
s_total = s_body + s_rotorHubs + s_blades;

% Signal used for time-frequency analysis
% Static point scatterers produce a strong zero-Doppler component.
if removeStaticMeanForTF
    s_tf = s_total - mean(s_total);
else
    s_tf = s_total;
end

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

        % Radial velocity under far-field approximation
        Rdot_tip = u_LOS.' * v_tip;

        % Monostatic Doppler frequency
        fd_tip(idx, :) = -(2/lambda) * Rdot_tip;
    end
end

%% 11. Plot : time-domain basband signal
figure('Name','Received Baseband Signal');
plot(t, real(s_total), 'LineWidth', 1.0); grid on;
xlabel('Time [s]');
ylabel('Real\{s_{total}(t)\}');
title('Total Baseband Return: Body + Rotor Hubs + Line-Integral Blades');

figure('Name','Blade-Only Baseband Signal');
plot(t, real(s_blades), 'LineWidth', 1.0); grid on;
xlabel('Time [s]');
ylabel('Real\{s_{blades}(t)\}');
title('Blade-Only Baseband Return');

%% 12. Plot: micro-Doppler spectrogram
winLen = 4096;
win = hamming(winLen, 'periodic');
noverlap = round(0.90 * winLen);
nfft = 32768;

[S, F, Tspec] = spectrogram(s_tf, win, noverlap, nfft, fs_slow, 'centered');

SdB = 20*log10(abs(S) ./ (max(abs(S(:))) + eps) + eps);

figure('Name','Micro-Doppler Spectrogram');
imagesc(Tspec, F, SdB);
axis xy;
ylim(1.5 * [-fd_max_approx fd_max_approx]);
clim([-40 0]);
grid on;
colorbar;
colormap(jet);
xlabel('Time [s]');
ylabel('Doppler frequency [Hz]');
title('Hovering Hexacopter Micro-Doppler Signature');

% % Optional analytical blade-tip guide curves
% hold on;
% plot(t, fd_tip, 'w--', 'LineWidth', 0.4);
% hold off;

%% 13. Plot: blade-only spectrogram
% [Sb, Fb, Tb] = spectrogram(s_blades, win, noverlap, nfft, fs_slow, 'centered');
% 
% Sb_dB = 20*log10(abs(Sb) ./ (max(abs(Sb(:))) + eps) + eps);
% 
% figure('Name','Blade-Only Micro-Doppler Spectrogram');
% imagesc(Tb, Fb, Sb_dB);
% axis xy;
% ylim(1.5 * [-fd_max_approx fd_max_approx]);
% clim([-40 0]);
% grid on;
% colorbar;
% colormap(jet);
% xlabel('Time [s]');
% ylabel('Doppler frequency [Hz]');
% title('Blade-Only Hovering Hexacopter Micro-Doppler Signature');

% hold on;
% plot(t, fd_tip, 'w--', 'LineWidth', 0.4);
% hold off;

%% 14. Plot: power spectrum
NfftSpec = 2^nextpow2(Nt);
Sfreq = fftshift(fft(s_tf, NfftSpec));
faxis = (-NfftSpec/2:NfftSpec/2-1) * fs_slow/NfftSpec;

figure('Name','Power Spectrum');
plot(faxis, 20*log10(abs(Sfreq)./(max(abs(Sfreq)) + eps) + eps), ...
    'LineWidth', 1.0);
grid on;
xlim(1.5 * [-fd_max_approx fd_max_approx]);
ylim([-80 5]);
xlabel('Doppler frequency [Hz]');
ylabel('Normalized magnitude [dB]');
title('Power Spectrum of Hovering Hexacopter Return');

%% 15. Plot: per-rotor blade returns
% figure('Name','Per-Rotor Blade Return Magnitudes');
% hold on; grid on;
% 
% for j = 1:NR
%     plot(t, abs(s_rotorBlades(j, :)), 'LineWidth', 0.8);
% end
% 
% xlabel('Time [s]');
% ylabel('|s_{rotor,j}(t)|');
% title('Per-Rotor Blade Return Magnitudes');
% legend(arrayfun(@(j) sprintf('Rotor %d', j), 1:NR, 'UniformOutput', false), ...
%     'Location', 'best');

%% 16. Plot: 3-D geometry
figure('Name','3-D Geometry');
hold on; grid on; axis equal;

% Radar Q
plot3(0, 0, 0, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
text(0, 0, 0, '  Q: Radar', 'FontSize', 10);

% Radar coordinate axes at Q
radarAxisLen = 5;
quiver3(0, 0, 0, radarAxisLen, 0, 0, ...
    'LineWidth', 1.0, 'MaxHeadSize', 1.0);
text(radarAxisLen, 0, 0, ' X');

quiver3(0, 0, 0, 0, radarAxisLen, 0, ...
    'LineWidth', 1.0, 'MaxHeadSize', 1.0);
text(0, radarAxisLen, 0, ' Y');

quiver3(0, 0, 0, 0, 0, radarAxisLen, ...
    'LineWidth', 1.0, 'MaxHeadSize', 1.0);
text(0, 0, radarAxisLen, ' Z');

% Hexacopter center O
plot3(O_radar(1), O_radar(2), O_radar(3), 'ro', ...
    'MarkerSize', 9, 'MarkerFaceColor', 'r');
text(O_radar(1), O_radar(2), O_radar(3), ...
    '  O: Hexacopter center', 'FontSize', 10);

% LOS from Q to O
plot3([0 O_radar(1)], [0 O_radar(2)], [0 O_radar(3)], ...
    'k--', 'LineWidth', 1.0);

% Reference/local coordinate axes at O
axisLen = 0.45;

quiver3(O_radar(1), O_radar(2), O_radar(3), axisLen, 0, 0, ...
    'LineWidth', 1.3, 'MaxHeadSize', 1.0);
text(O_radar(1)+axisLen, O_radar(2), O_radar(3), ' x / X''');

quiver3(O_radar(1), O_radar(2), O_radar(3), 0, axisLen, 0, ...
    'LineWidth', 1.3, 'MaxHeadSize', 1.0);
text(O_radar(1), O_radar(2)+axisLen, O_radar(3), ' y / Y''');

quiver3(O_radar(1), O_radar(2), O_radar(3), 0, 0, axisLen, ...
    'LineWidth', 1.3, 'MaxHeadSize', 1.0);
text(O_radar(1), O_radar(2), O_radar(3)+axisLen, ' z / Z''');

% Draw rotor hubs, arms, blade disks, and blades at t = 0
theta_circle = linspace(0, 2*pi, 200);

for j = 1:NR

    hub = O_radar + rotor_local(:, j);

    % Arm from O to rotor hub
    plot3([O_radar(1) hub(1)], ...
          [O_radar(2) hub(2)], ...
          [O_radar(3) hub(3)], ...
          'k-', 'LineWidth', 1.0);

    % Rotor hub
    plot3(hub(1), hub(2), hub(3), 'bo', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'b');

    if spin_sign(j) > 0
        dirText = 'CCW';
    else
        dirText = 'CW';
    end

    text(hub(1), hub(2), hub(3), ...
        sprintf('  R%d (%s)', j, dirText), 'FontSize', 9);

    % Rotor disk circle
    circle = hub + blade_tip * ...
        [cos(theta_circle);
         sin(theta_circle);
         zeros(size(theta_circle))];

    plot3(circle(1, :), circle(2, :), circle(3, :), ...
        ':', 'LineWidth', 0.8);

    % Rotation direction arrow on rotor disk
    arrow_theta = psi_rotor(j) + pi/4;
    p_arrow = hub + blade_tip * ...
        [cos(arrow_theta);
         sin(arrow_theta);
         0];

    tangent = spin_sign(j) * ...
        [-sin(arrow_theta);
          cos(arrow_theta);
          0];

    quiver3(p_arrow(1), p_arrow(2), p_arrow(3), ...
        0.03*tangent(1), 0.03*tangent(2), 0.03*tangent(3), ...
        'LineWidth', 1.2, 'MaxHeadSize', 2.0);

    % Blades at t = 0
    for k = 1:NB
        theta0 = phi0(j) + 2*pi*(k-1)/NB;

        u0 = [cos(theta0);
              sin(theta0);
              0];

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
title('Radar, Reference/Local Coordinates, Hexacopter Hubs, and Blades at t = 0');
view(35, 25);

%% 17. Plot: magnified top view around O
figure('Name','Magnified Top View Around Hexacopter Center');
hold on; grid on; axis equal;

% Local coordinate top-view, centered at O
plot(0, 0, 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
text(0, 0, '  O', 'FontSize', 10);

% Local x-y axes
quiver(0, 0, 0.5, 0, 'LineWidth', 1.2, 'MaxHeadSize', 1.0);
text(0.52, 0, 'x / X''');

quiver(0, 0, 0, 0.5, 'LineWidth', 1.2, 'MaxHeadSize', 1.0);
text(0, 0.52, 'y / Y''');

for j = 1:NR

    hub_local = rotor_local(:, j);

    % Arm
    plot([0 hub_local(1)], [0 hub_local(2)], ...
        'k-', 'LineWidth', 1.0);

    % Hub
    plot(hub_local(1), hub_local(2), 'bo', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'b');

    if spin_sign(j) > 0
        dirText = 'CCW';
    else
        dirText = 'CW';
    end

    text(hub_local(1), hub_local(2), ...
        sprintf('  R%d (%s)', j, dirText), 'FontSize', 9);

    % Rotor disk
    plot(hub_local(1) + blade_tip*cos(theta_circle), ...
         hub_local(2) + blade_tip*sin(theta_circle), ...
         ':', 'LineWidth', 0.8);

    % Blades at t = 0
    for k = 1:NB
        theta0 = phi0(j) + 2*pi*(k-1)/NB;

        u0 = [cos(theta0);
              sin(theta0)];

        p_root = hub_local(1:2) + blade_root * u0;
        p_tip  = hub_local(1:2) + blade_tip  * u0;

        plot([p_root(1) p_tip(1)], ...
             [p_root(2) p_tip(2)], ...
             'LineWidth', 2.0);
    end
end

xlabel('x / X'' [m]');
ylabel('y / Y'' [m]');
title('Magnified Top View: Hovering Hexacopter Geometry at t = 0');

%% 18. Summary print
fprintf('\n===== Hexacopter Simulation Summary =====\n');
fprintf('Rotor configuration: regular hexagon\n');
fprintf('NR = %d rotors\n', NR);
fprintf('NB = %d blades per rotor\n', NB);
fprintf('Rotor radius from O = %.3f m\n', rotor_radius);
fprintf('Blade root = %.3f m, tip = %.3f m, length = %.3f m\n', ...
    blade_root, blade_tip, blade_len);
fprintf('Hovering rotor rate = %.1f rev/s\n', Omega_rps);
fprintf('Approx. fd max = %.1f Hz\n', fd_max_approx);
fprintf('Blade-passage period = %.6f s\n', 1/(NB*Omega_rps));
fprintf('Slow-time fs = %.1f Hz, Tobs = %.3f s\n', fs_slow, Tobs);
fprintf('Spin signs:\n');
disp(spin_sign);