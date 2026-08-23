%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Geometry of vibrating point scatterer
%
% Point Scatterer:  Point scatterer location in Component coordinates:
%                   X = 0 m, Y = 0 m, Z = 0 m
%
% By H.J. Park
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

%% Parameters

T = 5;                          % Total time duration
Fs = 1e3;                       % Simulation sample rate
nt = T*Fs;                      % Total number of samples
dt = 1/Fs;                      
t = 0:dt:T-dt;                  % Time span
F0 = 5e9;                       % Center frequency
c = 299792458;
lambda = c/F0;                  % wavelength

Fv = 5;                         % Vibrating frequency
Av = 0.05;                      % Vibrating displacement
Phiv = 0;                       % Vibrating inital phase
Dirv = [1; 0; 0];               % Vibrating direction in Local coordinates
Dirv = Dirv/norm(Dirv);
R0 = [1; 0; 0];                 % initial position of target
V = [0; 0; 0];                  % Bulk velocity of Target
eulerInitDeg = [0; 0; 0];       % initial rotation of target
eulerRateDeg = [0; 0; 0];       % Bulk rotation of target
rcs = 0.8;                      % rcs of point scatterer, ideal

runGeometryAnimation = false;    % Enable 3-D geometry animation
animationFrameRate = 30;        % Display frame rate [frames/s]
animationPlaybackSpeed = 1;     % 1: real time, 2: twice as fast


%% Geometry

% Radar position in Radar coordinates
RadarPos = [0; 0; 0];                       

% Translation of target center about radar
RefRadarTrans = [R0(1)+V(1)*t; ...
                 R0(2)+V(2)*t; ...
                 R0(3)+V(3)*t];

% Rotation of target about radar
yawchange = eulerRateDeg(1)*t;
pitchchange = eulerRateDeg(2)*t;
rollchange = eulerRateDeg(3)*t;
yaw = eulerInitDeg(1)+yawchange;
pitch = eulerInitDeg(2)+pitchchange;
roll = eulerInitDeg(3)+rollchange;
N = numel(t);
LocRefRot = zeros(3,3,N);
for k=1:N
    cy = cosd(yaw(k)); sy = sind(yaw(k));
    cp = cosd(pitch(k)); sp = sind(pitch(k));
    cr = cosd(roll(k)); sr = sind(roll(k));
    Rz = [cy, -sy, 0; sy, cy, 0; 0, 0, 1];
    Ry = [cp, 0, sp; 0, 1, 0; -sp, 0, cp];
    Rx = [1, 0, 0; 0, cr, -sr; 0, sr, cr];
    LocRefRot(:,:,k) = Rz * Ry * Rx;
end

% Translation of component about target center
vx = Dirv(1)*Av*sin(2*pi*Fv*t+Phiv);
vy = Dirv(2)*Av*sin(2*pi*Fv*t+Phiv);
vz = Dirv(3)*Av*sin(2*pi*Fv*t+Phiv);
CompLocTrans = [0+vx; ...
                0+vy; ...
                0+vz];                      

% Rotation of component about target
vyaw = 0;
vpitch = 0;
vroll = 0;
vcy = cosd(vyaw); vsy = sind(vyaw);
vcp = cosd(vpitch); vsp = sind(vpitch);
vcr = cosd(vroll); vsr = sind(vroll);
vRz = [vcy, -vsy, 0; vsy, vcy, 0; 0, 0, 1];
vRy = [vcp, 0, vsp; 0, 1, 0; -vsp, 0, vcp];
vRx = [1, 0, 0; 0, vcr, -vsr; 0, vsr, vcr];
CompLocRot = vRz * vRy * vRx;               

% Position of point scatterer in Component coordinates
PointCompPos = [0; 0; 0];

%% Simulation Scenario

data = zeros(1,N);
PointPos = zeros(3,N);

for k=1:N
    % position of point scatterer in Radar coordinates
    PointPos(:,k) = RefRadarTrans(:,k) + LocRefRot(:,:,k) * ...
                (CompLocTrans(:,k) + CompLocRot * PointCompPos);
    % distance from radar to point scatterer
    r_dist(:,k) = PointPos(:,k) - RadarPos(:);
    distance(k) = sqrt(r_dist(1,k).^2+r_dist(2,k).^2+r_dist(3,k).^2);
    amp(k) = sqrt(rcs);
    PHs = amp(k)*exp(-1j*4*pi*distance(k)/lambda);
    data(k) = PHs;
end

%% Plot results

% Micro-Doppler Signature
x = data;
figure(1)
colormap(jet)
spectrogram(x,kaiser(64,10),60,512,Fs,'centered','yaxis');
clim = get(gca,'CLim');
set(gca,'CLim',clim(2) + [-50 0]);
drawnow

%% Geometry animation

if runGeometryAnimation
    % Only a subset of the simulation samples is rendered. The signal and
    % range calculations above still use every sample.
    animationStep = max(1, round(Fs*animationPlaybackSpeed/animationFrameRate));
    animationIndex = 1:animationStep:N;
    if animationIndex(end) ~= N
        animationIndex(end+1) = N;
    end

    % Full-scene limits include the radar, target center, and scatterer.
    geometryPoints = [RadarPos, RefRadarTrans, PointPos];
    geometryMin = min(geometryPoints, [], 2);
    geometryMax = max(geometryPoints, [], 2);
    geometrySpan = geometryMax - geometryMin;
    sceneScale = max(geometrySpan);
    if sceneScale == 0
        sceneScale = 1;
    end
    sceneMargin = 0.10*sceneScale;
    targetAxisLength = 0.08*sceneScale;
    radarAxisLength = targetAxisLength;

    % The zoomed view removes target-center translation so that the small
    % vibration remains visible even when the radar range is much larger.
    PointRelTarget = PointPos - RefRadarTrans;
    maxRelativeExtent = max(abs(PointRelTarget), [], 'all');
    zoomRadius = max([2.0*maxRelativeExtent, 2.0*Av, 0.05]);
    zoomAxisLength = 0.65*zoomRadius;

    geometryFigure = figure('Name', 'Vibrating Point Scatterer Geometry', ...
        'Color', 'w', 'Position', [100 100 1450 650]);
    animationLayout = tiledlayout(geometryFigure, 1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    %% Full radar-target geometry
    sceneAxes = nexttile(animationLayout, 1);
    hold(sceneAxes, 'on');
    grid(sceneAxes, 'on');
    axis(sceneAxes, 'equal');
    view(sceneAxes, 35, 25);
    xlabel(sceneAxes, 'Radar X [m]');
    ylabel(sceneAxes, 'Radar Y [m]');
    zlabel(sceneAxes, 'Radar Z [m]');
    xlim(sceneAxes, [geometryMin(1)-sceneMargin, geometryMax(1)+sceneMargin]);
    ylim(sceneAxes, [geometryMin(2)-sceneMargin, geometryMax(2)+sceneMargin]);
    zlim(sceneAxes, [geometryMin(3)-sceneMargin, geometryMax(3)+sceneMargin]);

    radarHandle = plot3(sceneAxes, RadarPos(1), RadarPos(2), RadarPos(3), ...
        'kp', 'MarkerFaceColor', [1.0 0.8 0.1], 'MarkerSize', 13, ...
        'DisplayName', 'Radar');
    centerPathHandle = plot3(sceneAxes, RefRadarTrans(1,:), ...
        RefRadarTrans(2,:), RefRadarTrans(3,:), 'k--', 'LineWidth', 1.0, ...
        'DisplayName', 'Target-center path');
    centerHandle = plot3(sceneAxes, RefRadarTrans(1,1), ...
        RefRadarTrans(2,1), RefRadarTrans(3,1), 'ko', ...
        'MarkerFaceColor', [0.2 0.7 0.9], 'MarkerSize', 8, ...
        'DisplayName', 'Target center');
    pointHandle = plot3(sceneAxes, PointPos(1,1), PointPos(2,1), ...
        PointPos(3,1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7, ...
        'DisplayName', 'Point scatterer');
    pointTraceHandle = animatedline(sceneAxes, 'Color', [0.85 0.2 0.2], ...
        'LineWidth', 1.2, 'DisplayName', 'Scatterer path');
    losHandle = plot3(sceneAxes, [RadarPos(1), PointPos(1,1)], ...
        [RadarPos(2), PointPos(2,1)], [RadarPos(3), PointPos(3,1)], ...
        'm-', 'LineWidth', 1.0, 'DisplayName', 'Radar LOS');
    offsetHandle = plot3(sceneAxes, [RefRadarTrans(1,1), PointPos(1,1)], ...
        [RefRadarTrans(2,1), PointPos(2,1)], ...
        [RefRadarTrans(3,1), PointPos(3,1)], 'c-', 'LineWidth', 2.0, ...
        'DisplayName', 'Center-to-scatterer');

    % Static radar-coordinate axes.
    quiver3(sceneAxes, RadarPos(1), RadarPos(2), RadarPos(3), ...
        radarAxisLength, 0, 0, 0, 'Color', [0.7 0 0], ...
        'LineStyle', ':', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    quiver3(sceneAxes, RadarPos(1), RadarPos(2), RadarPos(3), ...
        0, radarAxisLength, 0, 0, 'Color', [0 0.55 0], ...
        'LineStyle', ':', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    quiver3(sceneAxes, RadarPos(1), RadarPos(2), RadarPos(3), ...
        0, 0, radarAxisLength, 0, 'Color', [0 0.25 0.8], ...
        'LineStyle', ':', 'LineWidth', 1.3, 'HandleVisibility', 'off');
    text(sceneAxes, RadarPos(1)+radarAxisLength, RadarPos(2), RadarPos(3), ...
        'X_R', 'Color', [0.7 0 0]);
    text(sceneAxes, RadarPos(1), RadarPos(2)+radarAxisLength, RadarPos(3), ...
        'Y_R', 'Color', [0 0.55 0]);
    text(sceneAxes, RadarPos(1), RadarPos(2), RadarPos(3)+radarAxisLength, ...
        'Z_R', 'Color', [0 0.25 0.8]);

    initialCenter = RefRadarTrans(:,1);
    initialRotation = LocRefRot(:,:,1);
    targetXAxis = quiver3(sceneAxes, initialCenter(1), initialCenter(2), ...
        initialCenter(3), targetAxisLength*initialRotation(1,1), ...
        targetAxisLength*initialRotation(2,1), ...
        targetAxisLength*initialRotation(3,1), 0, 'r', 'LineWidth', 2.0, ...
        'DisplayName', 'Target X_L');
    targetYAxis = quiver3(sceneAxes, initialCenter(1), initialCenter(2), ...
        initialCenter(3), targetAxisLength*initialRotation(1,2), ...
        targetAxisLength*initialRotation(2,2), ...
        targetAxisLength*initialRotation(3,2), 0, 'g', 'LineWidth', 2.0, ...
        'DisplayName', 'Target Y_L');
    targetZAxis = quiver3(sceneAxes, initialCenter(1), initialCenter(2), ...
        initialCenter(3), targetAxisLength*initialRotation(1,3), ...
        targetAxisLength*initialRotation(2,3), ...
        targetAxisLength*initialRotation(3,3), 0, 'b', 'LineWidth', 2.0, ...
        'DisplayName', 'Target Z_L');

    legend(sceneAxes, [radarHandle, centerPathHandle, centerHandle, ...
        pointHandle, pointTraceHandle, losHandle, offsetHandle, ...
        targetXAxis, targetYAxis, targetZAxis], 'Location', 'best');

    %% Target-centered zoomed geometry
    zoomAxes = nexttile(animationLayout, 2);
    hold(zoomAxes, 'on');
    grid(zoomAxes, 'on');
    axis(zoomAxes, 'equal');
    view(zoomAxes, 35, 25);
    xlabel(zoomAxes, 'Relative X_R [m]');
    ylabel(zoomAxes, 'Relative Y_R [m]');
    zlabel(zoomAxes, 'Relative Z_R [m]');
    xlim(zoomAxes, zoomRadius*[-1 1]);
    ylim(zoomAxes, zoomRadius*[-1 1]);
    zlim(zoomAxes, zoomRadius*[-1 1]);
    title(zoomAxes, 'Target-centered vibration (translation removed)');

    zoomCenterHandle = plot3(zoomAxes, 0, 0, 0, 'ko', ...
        'MarkerFaceColor', [0.2 0.7 0.9], 'MarkerSize', 8, ...
        'DisplayName', 'Target center');
    zoomPointHandle = plot3(zoomAxes, PointRelTarget(1,1), ...
        PointRelTarget(2,1), PointRelTarget(3,1), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
        'DisplayName', 'Point scatterer');
    zoomTraceHandle = animatedline(zoomAxes, 'Color', [0.85 0.2 0.2], ...
        'LineWidth', 1.5, 'DisplayName', 'Relative trajectory');
    zoomOffsetHandle = plot3(zoomAxes, [0, PointRelTarget(1,1)], ...
        [0, PointRelTarget(2,1)], [0, PointRelTarget(3,1)], ...
        'c-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    zoomXAxis = quiver3(zoomAxes, 0, 0, 0, ...
        zoomAxisLength*initialRotation(1,1), ...
        zoomAxisLength*initialRotation(2,1), ...
        zoomAxisLength*initialRotation(3,1), 0, 'r', 'LineWidth', 2.0, ...
        'DisplayName', 'Target X_L');
    zoomYAxis = quiver3(zoomAxes, 0, 0, 0, ...
        zoomAxisLength*initialRotation(1,2), ...
        zoomAxisLength*initialRotation(2,2), ...
        zoomAxisLength*initialRotation(3,2), 0, 'g', 'LineWidth', 2.0, ...
        'DisplayName', 'Target Y_L');
    zoomZAxis = quiver3(zoomAxes, 0, 0, 0, ...
        zoomAxisLength*initialRotation(1,3), ...
        zoomAxisLength*initialRotation(2,3), ...
        zoomAxisLength*initialRotation(3,3), 0, 'b', 'LineWidth', 2.0, ...
        'DisplayName', 'Target Z_L');

    legend(zoomAxes, [zoomCenterHandle, zoomPointHandle, zoomTraceHandle, ...
        zoomXAxis, zoomYAxis, zoomZAxis], 'Location', 'best');

    %% Update animated geometry
    for frameNumber = 1:numel(animationIndex)
        if ~isgraphics(geometryFigure)
            break
        end

        k = animationIndex(frameNumber);
        centerPosition = RefRadarTrans(:,k);
        pointPosition = PointPos(:,k);
        relativePosition = PointRelTarget(:,k);
        targetRotation = LocRefRot(:,:,k);

        set(centerHandle, 'XData', centerPosition(1), ...
            'YData', centerPosition(2), 'ZData', centerPosition(3));
        set(pointHandle, 'XData', pointPosition(1), ...
            'YData', pointPosition(2), 'ZData', pointPosition(3));
        set(losHandle, 'XData', [RadarPos(1), pointPosition(1)], ...
            'YData', [RadarPos(2), pointPosition(2)], ...
            'ZData', [RadarPos(3), pointPosition(3)]);
        set(offsetHandle, 'XData', [centerPosition(1), pointPosition(1)], ...
            'YData', [centerPosition(2), pointPosition(2)], ...
            'ZData', [centerPosition(3), pointPosition(3)]);
        addpoints(pointTraceHandle, pointPosition(1), ...
            pointPosition(2), pointPosition(3));

        set(targetXAxis, 'XData', centerPosition(1), ...
            'YData', centerPosition(2), 'ZData', centerPosition(3), ...
            'UData', targetAxisLength*targetRotation(1,1), ...
            'VData', targetAxisLength*targetRotation(2,1), ...
            'WData', targetAxisLength*targetRotation(3,1));
        set(targetYAxis, 'XData', centerPosition(1), ...
            'YData', centerPosition(2), 'ZData', centerPosition(3), ...
            'UData', targetAxisLength*targetRotation(1,2), ...
            'VData', targetAxisLength*targetRotation(2,2), ...
            'WData', targetAxisLength*targetRotation(3,2));
        set(targetZAxis, 'XData', centerPosition(1), ...
            'YData', centerPosition(2), 'ZData', centerPosition(3), ...
            'UData', targetAxisLength*targetRotation(1,3), ...
            'VData', targetAxisLength*targetRotation(2,3), ...
            'WData', targetAxisLength*targetRotation(3,3));

        set(zoomPointHandle, 'XData', relativePosition(1), ...
            'YData', relativePosition(2), 'ZData', relativePosition(3));
        set(zoomOffsetHandle, 'XData', [0, relativePosition(1)], ...
            'YData', [0, relativePosition(2)], ...
            'ZData', [0, relativePosition(3)]);
        addpoints(zoomTraceHandle, relativePosition(1), ...
            relativePosition(2), relativePosition(3));

        set(zoomXAxis, 'UData', zoomAxisLength*targetRotation(1,1), ...
            'VData', zoomAxisLength*targetRotation(2,1), ...
            'WData', zoomAxisLength*targetRotation(3,1));
        set(zoomYAxis, 'UData', zoomAxisLength*targetRotation(1,2), ...
            'VData', zoomAxisLength*targetRotation(2,2), ...
            'WData', zoomAxisLength*targetRotation(3,2));
        set(zoomZAxis, 'UData', zoomAxisLength*targetRotation(1,3), ...
            'VData', zoomAxisLength*targetRotation(2,3), ...
            'WData', zoomAxisLength*targetRotation(3,3));

        title(sceneAxes, sprintf( ...
            'Full geometry: t = %.3f s, range = %.3f m', ...
            t(k), distance(k)));
        drawnow;
        pause(1/animationFrameRate);
    end
end
