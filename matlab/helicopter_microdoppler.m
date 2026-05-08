%% Helicopter micro-Doppler simulation in 3D
% based on Chen, "The micro-Doppler effect in radar"
%
% Models
%   - Blades    : line-integral model
%   - Body/Hub  : point-scattering model
clear; clc; close all;

%% Radar Parameters
c = 3e8;                % speed of light
lambda = 0.06;          % wavelength (C-band)
fc = c / lambda;        % carrier frequency

%% Helicopter Geometry
R0 = 700;               % radar-to-helicopter range
az = deg2rad(20);       % azimuth angle
el = deg2rad(45);       % elevation angle

%% Rotor Parameters
Nb = 3;                 % number of blades per rotor
L = 6.5;                % blade length
frot = 2;               % rotor rotation rate
Omega = 2*pi*frot;      % angular speed

%% Scatterer Amplitudes
A_body = 1.0;           % body RCS amplitude
A_hub = 0.5;            % hub RCS amplitude

%% Sampling
v_tip = Omega * L;
fd_max = 2*v_tip*cos(el)/lambda;

fprintf('Tip speed             : %.2f m/s\n', v_tip);
fprintf('Maximum micro-Doppler : %.2f Hz\n', fd_max);
fprintf('Blade-flash period    : %.2f ms\n', 1e3/(Nb*frot));

fs = 8e3;               % sampling frequency
T = 1.0;                % observation interval
t = (0:1/fs:T-1/fs);
N = length(t);

%% 3D Geometry
% Radar at the origin. LOS unit vector points from radar to helicopter.
u_LOS = [cos(el)*cos(az); cos(el)*sin(az); sin(el)];

% Rotor center, hub, body
P_body = R0 * u_LOS;
P_hub = P_body + [0; 0; 1];         % hub 1m above the body
P_rotor = P_hub;                    % rotor centre coincides with hub

%% Signal Synthesis

% body - point scatterer
R_body = norm(P_body);
s_body = A_body * exp(-1j*4*pi*R_body/lambda) * ones(1,N);

% hub - point scatterer
R_hub = norm(P_hub);
s_hub = A_hub * exp(-1j*4*pi*R_hub/lambda) * ones(1,N);

% blades - line-integral model
s_blades = zeros(1,N);
for k = 0:Nb-1
    phi_k = 2*pi*k/Nb;                      % initial blade phase
    theta_k = Omega*t + phi_k;              % time-varying angle
    proj = cos(el) .* cos(theta_k - az);    % LOS projection, far-field approximation
    arg = (2*L/lambda) .* proj;             % sinc argument
    mag = L * sinc(arg);                    % MATLAB sinc: sin(pi x)/(pi x)
    phase = exp(-1j*4*pi/lambda * (R0 + sin(el) + (L/2).*proj));
    s_blades = s_blades + mag .* phase;
end

% total received complex-baseband signal
s = s_body + s_hub + s_blades;

%% Plot 3D geometry
figure('Name','3D Geometry','Color','w'); hold on; grid on; axis equal;
 
% Radar position (at origin)
hRadar = plot3(0,0,0,'r^','MarkerSize',12,'MarkerFaceColor','r');
text(0,0,0,'  Radar','FontWeight','bold');
 
% Line of sight (radar -> rotor)
plot3([0 P_body(1)],[0 P_body(2)],[0 P_body(3)],'r--','HandleVisibility','off');
 
% Rotor blades (line-integral model) at t = 0
for k = 0:Nb-1
    phi_k = 2*pi*k/Nb;
    bt    = P_rotor + L*[cos(phi_k); sin(phi_k); 0];
    if k == 0
        hBlade = plot3([P_rotor(1) bt(1)],[P_rotor(2) bt(2)],[P_rotor(3) bt(3)], ...
                       'b-','LineWidth',2);
    else
        plot3([P_rotor(1) bt(1)],[P_rotor(2) bt(2)],[P_rotor(3) bt(3)], ...
              'b-','LineWidth',2,'HandleVisibility','off');
    end
end
 
% Mast (line connecting hub to body)
plot3([P_hub(1) P_body(1)],[P_hub(2) P_body(2)],[P_hub(3) P_body(3)], ...
      'k-','LineWidth',1.5,'HandleVisibility','off');
 
% Rotor hub (point scatterer)
hHub  = plot3(P_hub(1), P_hub(2), P_hub(3), ...
              'go','MarkerSize',10,'MarkerFaceColor','g','LineWidth',1.2);
text(P_hub(1),P_hub(2),P_hub(3),'  Hub','FontWeight','bold');
 
% Helicopter body (point scatterer)
hBody = plot3(P_body(1),P_body(2),P_body(3), ...
              'ks','MarkerSize',12,'MarkerFaceColor',[0.3 0.3 0.3]);
text(P_body(1),P_body(2),P_body(3),'  Body','FontWeight','bold');
 
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('3D Geometry: Radar, Body, Hub, and Rotor Blades (t = 0)');
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