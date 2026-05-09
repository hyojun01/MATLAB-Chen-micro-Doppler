%% Helicopter micro-Doppler Simulation in 3D
%
% Models
%   body / hub : point-scattering model
%   blades : line-integral model
%
% 
clear; clc; close all;

%% radar parameters
c = 3e8;                % speed of light
fc = 5e9;               % carrier frequency
lambda = c / fc;        % wavelength

%% helicopter geometry & motion
R0 = 700;                                                       % initial range to helicopter
az = deg2rad(20);                                               % initial azimuth angle
el = deg2rad(45);                                               % initial elevation angle
u_LOS_0 = [cos(el)*cos(az); cos(el)*sin(az); sin(el)];          % initial unit vector from radar to body
P_body_0 = R0 * u_LOS_0;                                        % initial position of body in radar coordinates
v = [50; 0; 0];                                                 % translation velocity

% initial body rotation
theta_tilt_deg = 7;
Rx = @(a) [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
Ry = @(a) [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
R_init = Rz(90) * Rx(theta_tilt_deg) * Rz(-90);

%% rotor parameters
Nb = 2;                 % number of blades
L = 5.03;               % blade length
frot = 6.8;             % rotor frequency
Omega = 2*pi*frot;      % rotor angular speed

% body position in local coordinates
r_body = [0; 0; 0];

% rotor hub position in local coordinates
r_hub = [0; 0; 1.0];

% initial blade phases
phi_k = 2*pi*(0:Nb-1).'/Nb;

%% scatterer amplitudes
A_body = 1.0;           % body RCS amplitude
A_hub = 0.5;            % rotor hub RCS amplitude
A_blade = 1.0;          % rotor blade RCS amplitude

%% Sampling
v_tip = Omega * L;
fd_max = 2 * v_tip / lambda;        % max blade-tip Doppler

fprintf('Tip speed             : %.2f m/s\n', v_tip);
fprintf('Maximum micro-Doppler : %.2f Hz\n', fd_max);
fprintf('Blade-flash period    : %.2f ms\n', 1e3/(Nb*frot));

fs = 25e3;              % PRF / sampling frequency
T = 1.0;                % observation interval
t = 0:1/fs:T-1/fs;
N = length(t);

%% 3D trajectory and line of sight
P_body = P_body_0 + v * t;              % body trajectory
R_body = sqrt(sum(P_body.^2, 1));       % range from radar to body
u_LOS = P_body / R_body;                % unit vector from radar to body

%% signal synthesis

% body - point scatterer
s_body = A_body * exp(-1j*4*pi/lambda * R_body);

% hub - point scatterer
r_hub_ref = R_init * r_hub;             % rotor hub position in reference coordinates
A = (u_LOS.' * r_hub_ref).';
R_hub = R_body + A;
s_hub = A_hub * exp(-1j*4*pi/lambda * R_hub);

% blades - line-integral model
s_blades = zeros(1,N);
for k = 1:Nb
    theta_k = Omega*t + phi_k(k);
    b_local = [cos(theta_k); sin(theta_k); zeros(1,N)];
    b_ref = R_init * b_local;
    proj = sum(u_LOS .* b_ref, 1);

    arg = (2*L/lambda) .* proj;         % sinc argument
    mag = L * sinc(arg);
    phase = exp(-1j*4*pi/lambda * (R_body + A + (L/2).*proj));
    s_blades = s_blades + A_blade * mag .* phase;
end

% total received complex-baseband signal
s = s_body + s_hub + s_blades;

%% Plot 3D geometry 
P_body0 = P_body(:,1);
P_hub0  = P_body0 + R_init * r_hub;
 
figure('Name','3D Geometry','Color','w'); hold on; grid on; axis equal;
 
% Radar (origin)
hRadar = plot3(0,0,0,'r^','MarkerSize',12,'MarkerFaceColor','r');
text(0,0,0,'  Radar','FontWeight','bold');
 
% Line of sight (radar -> body)
plot3([0 P_body0(1)], [0 P_body0(2)], [0 P_body0(3)], ...
      'r--','HandleVisibility','off');
 
% Velocity vector (scaled for visibility)
quiver3(P_body0(1),P_body0(2),P_body0(3), ...
        v(1)*0.4, v(2)*0.4, v(3)*0.4, ...
        0,'m','LineWidth',1.5,'MaxHeadSize',0.6,'HandleVisibility','off');
text(P_body0(1)+v(1)*0.4, P_body0(2), P_body0(3), '  v', 'Color','m');
 
% Rotor blades at t = 0  (R_init applied)
for k = 1:Nb
    bt_local = r_hub + L*[cos(phi_k(k)); sin(phi_k(k)); 0];
    bt_ref   = P_body0 + R_init * bt_local;
    if k == 1
        hBlade = plot3([P_hub0(1) bt_ref(1)],[P_hub0(2) bt_ref(2)], ...
                       [P_hub0(3) bt_ref(3)],'b-','LineWidth',2);
    else
        plot3([P_hub0(1) bt_ref(1)],[P_hub0(2) bt_ref(2)], ...
              [P_hub0(3) bt_ref(3)],'b-','LineWidth',2,'HandleVisibility','off');
    end
end
 
% Mast (line connecting hub to body)
plot3([P_hub0(1) P_body0(1)],[P_hub0(2) P_body0(2)],[P_hub0(3) P_body0(3)], ...
      'k-','LineWidth',1.5,'HandleVisibility','off');
 
% Rotor hub (point scatterer)
hHub  = plot3(P_hub0(1),P_hub0(2),P_hub0(3), ...
              'go','MarkerSize',10,'MarkerFaceColor','g','LineWidth',1.2);
text(P_hub0(1),P_hub0(2),P_hub0(3),'  Hub','FontWeight','bold');
 
% Helicopter body (point scatterer)
hBody = plot3(P_body0(1),P_body0(2),P_body0(3), ...
              'ks','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
text(P_body0(1),P_body0(2),P_body0(3),'  Body','FontWeight','bold');
 
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('3D Geometry at t = 0  (forward pitch %d°)', theta_tilt_deg));
legend([hRadar hBlade hHub hBody], ...
       {'Radar','Blade (line scatterer)','Hub (point scatterer)', ...
        'Body (point scatterer)'},'Location','best');
view(35,25);

%% Plot time domain
figure('Name','Received Signal','Color','w');
subplot(2,1,1); plot(t, real(s)); grid on; xlim([0 0.5]);
xlabel('Time [s]'); ylabel('Re\{s(t)\}');
title('Real part of complex-baseband signal');
subplot(2,1,2); plot(t, abs(s));  grid on; xlim([0 0.5]);
xlabel('Time [s]'); ylabel('|s(t)|');
title('Magnitude (note the periodic blade flashes)');

%% STFT / Spectrogram
win_len = 256;
nfft    = 2048;
overlap = win_len - 4;
[S,F,Ts] = spectrogram(s, hamming(win_len), overlap, nfft, fs, 'centered');
S_dB     = 20*log10(abs(S) + eps);
 
figure('Name','Micro-Doppler Spectrogram','Color','w');
imagesc(Ts, F, S_dB); axis xy;
colormap(jet); cb = colorbar; ylabel(cb,'Magnitude [dB]');
clim([max(S_dB(:))-50, max(S_dB(:))]);
xlabel('Time [s]'); ylabel('Doppler frequency [Hz]');
ylim([-1.1*fd_max, 1.1*fd_max]);
title('Helicopter Micro-Doppler Signature (STFT)');
 
% bulk-Doppler reference line (negative because target recedes)
v_radial0 = (v.') * u_LOS(:,1);
fd_body0  = -2*v_radial0/lambda;
hold on;
yline(fd_body0, 'w--', 'LineWidth', 1.4);