%% Vibration Micro-Doppler Baseband Signal
% Monostatic radar
%
% Positive displacement means that the target moves away from the radar.
%
% x(t)     = A_vib * sin(2*pi*f_vib*t)
% phi(t)   = -(4*pi/lambda) * x(t)
% sBB(t)   = exp(1j*phi(t))
% f_inst   = (1/(2*pi)) * d(phi)/dt
%          = -(2/lambda) * velocity(t)

clear;
close all;
clc;

%% 1. User parameters

centerFrequency   = 24e9;   % Radar center frequency [Hz]
vibrationRate     = 5;      % Vibration frequency [Hz]
vibrationAmplitude = 1e-2;  % Peak vibration amplitude [m]

numberOfPeriods   = 5;      % Number of vibration periods to display
runAnimation      = true;   % Animate phase motion on the complex plane

%% 2. Radar and vibration parameters

c = 299792458;                          % Speed of light [m/s]
lambda = c / centerFrequency;           % Wavelength [m]

omegaVib = 2*pi*vibrationRate;          % Vibration angular frequency [rad/s]

% Peak phase deviation
phaseDeviation = 4*pi*vibrationAmplitude/lambda;

% Maximum instantaneous micro-Doppler frequency
maxInstantaneousFrequency = ...
    4*pi*vibrationAmplitude*vibrationRate/lambda;

% Automatically select a sufficiently high simulation sampling frequency.
% This sampling frequency is only for numerical simulation and is not the
% RF ADC sampling frequency.
samplingFrequency = max(5e3, ...
    100*(vibrationRate + maxInstantaneousFrequency));
% samplingFrequency = 5e3;

simulationTime = numberOfPeriods / vibrationRate;

time = (0:1/samplingFrequency:simulationTime).';

%% 3. Vibration motion

displacement = vibrationAmplitude * sin(omegaVib*time);

velocity = vibrationAmplitude * omegaVib * cos(omegaVib*time);

%% 4. Baseband phase and complex signal

% For a monostatic radar, the round-trip phase is 4*pi*R/lambda.
% The constant phase associated with the mean target range is removed.
phaseExact = -(4*pi/lambda) * displacement;

basebandSignal = exp(1j*phaseExact);

% Wrapped phase: range [-pi, pi]
phaseWrapped = angle(basebandSignal);

% Recovered continuous phase
phaseUnwrapped = unwrap(phaseWrapped);

%% 5. Instantaneous frequency

% Analytical result
instantaneousFrequencyExact = -(2/lambda) * velocity;

% Equivalent closed-form expression
instantaneousFrequencyClosedForm = ...
    -maxInstantaneousFrequency * cos(omegaVib*time);

% Numerical estimate from the baseband phase
instantaneousFrequencyNumerical = ...
    gradient(phaseUnwrapped, 1/samplingFrequency)/(2*pi);

%% 6. Print calculated parameters

fprintf('====================================================\n');
fprintf('Vibration micro-Doppler parameters\n');
fprintf('====================================================\n');
fprintf('Center frequency              : %.3f GHz\n', ...
    centerFrequency/1e9);
fprintf('Wavelength                    : %.6f m\n', lambda);
fprintf('Vibration rate                : %.3f Hz\n', ...
    vibrationRate);
fprintf('Vibration amplitude           : %.3f mm\n', ...
    vibrationAmplitude*1e3);
fprintf('Peak-to-peak displacement     : %.3f mm\n', ...
    2*vibrationAmplitude*1e3);
fprintf('Peak phase deviation          : %.6f rad\n', ...
    phaseDeviation);
fprintf('Peak phase deviation          : %.3f deg\n', ...
    rad2deg(phaseDeviation));
fprintf('Maximum instantaneous freq.   : %.6f Hz\n', ...
    maxInstantaneousFrequency);
fprintf('Simulation sampling frequency : %.3f kHz\n', ...
    samplingFrequency/1e3);
fprintf('====================================================\n');

%% 7. Plot results

figure('Name', 'Vibration Micro-Doppler Baseband Signal', ...
    'Position', [100 100 1300 820]);

layout = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% 7-1. Complex baseband signal on the unit circle

nexttile;

circlePhase = linspace(0, 2*pi, 1000);

plot(cos(circlePhase), sin(circlePhase), '--', ...
    'LineWidth', 1);
hold on;

% Reduce the number of displayed points when necessary.
displayStep = max(1, floor(length(time)/1500));
displayIndex = 1:displayStep:length(time);

scatter(real(basebandSignal(displayIndex)), ...
        imag(basebandSignal(displayIndex)), ...
        18, time(displayIndex), 'filled');

plot(real(basebandSignal(1)), ...
     imag(basebandSignal(1)), ...
     'o', ...
     'MarkerSize', 9, ...
     'LineWidth', 2);

axis equal;
xlim([-1.15 1.15]);
ylim([-1.15 1.15]);
grid on;

xlabel('In-phase component');
ylabel('Quadrature component');
title('Baseband phase motion on the unit circle');

colorbarHandle = colorbar;
colorbarHandle.Label.String = 'Time [s]';

legend('Unit circle', ...
       'Time-varying baseband signal', ...
       'Starting point', ...
       'Location', 'best');

%% 7-2. Baseband signal in the time domain

nexttile;

plot(time, real(basebandSignal), ...
    'LineWidth', 1.3);
hold on;

plot(time, imag(basebandSignal), ...
    'LineWidth', 1.3);

grid on;
xlabel('Time [s]');
ylabel('Normalized amplitude');
title('Complex baseband signal versus time');

legend('Real\{s_{BB}(t)\}', ...
       'Imag\{s_{BB}(t)\}', ...
       'Location', 'best');

%% 7-3. Vibration displacement and phase

nexttile;

yyaxis left;

plot(time, displacement*1e3, ...
    'LineWidth', 1.3);

ylabel('Displacement [mm]');

yyaxis right;

plot(time, phaseExact, ...
    'LineWidth', 1.3);

ylabel('Baseband phase [rad]');

grid on;
xlabel('Time [s]');
title('Vibration displacement and baseband phase');

legend('Displacement x(t)', ...
       'Phase \phi(t)', ...
       'Location', 'best');

%% 7-4. Instantaneous frequency

nexttile;

plot(time, instantaneousFrequencyExact, ...
    'LineWidth', 1.5);
hold on;

plot(time, instantaneousFrequencyNumerical, '--', ...
    'LineWidth', 1.1);

yline(maxInstantaneousFrequency, ':', ...
    'Maximum');
yline(-maxInstantaneousFrequency, ':', ...
    'Minimum');

grid on;
xlabel('Time [s]');
ylabel('Instantaneous frequency [Hz]');
title('Vibration-induced instantaneous frequency');

legend('Analytical result', ...
       'Numerical phase derivative', ...
       'Location', 'best');

title(layout, sprintf([ ...
    'Vibration Micro-Doppler: f_c = %.3f GHz, ' ...
    'f_v = %.3f Hz, A_v = %.3f mm'], ...
    centerFrequency/1e9, ...
    vibrationRate, ...
    vibrationAmplitude*1e3));

%% 8. Validate analytical expressions

frequencyDifference = max(abs( ...
    instantaneousFrequencyExact - ...
    instantaneousFrequencyClosedForm));

fprintf('\nMaximum difference between analytical forms: %.3e Hz\n', ...
    frequencyDifference);

%% 9. Optional animation of phase motion on the complex plane

if runAnimation

    figure('Name', 'Animated Baseband Phase Motion', ...
        'Position', [250 150 750 700]);

    plot(cos(circlePhase), sin(circlePhase), '--', ...
        'LineWidth', 1);
    hold on;

    trajectory = animatedline( ...
        'LineWidth', 1.3, ...
        'MaximumNumPoints', 1000);

    movingPoint = plot(real(basebandSignal(1)), ...
                       imag(basebandSignal(1)), ...
                       'o', ...
                       'MarkerSize', 10, ...
                       'MarkerFaceColor', 'auto', ...
                       'LineWidth', 1.5);

    phaseVector = quiver(0, 0, ...
        real(basebandSignal(1)), ...
        imag(basebandSignal(1)), ...
        0, ...
        'LineWidth', 1.5);

    axis equal;
    xlim([-1.2 1.2]);
    ylim([-1.2 1.2]);
    grid on;

    xlabel('In-phase component');
    ylabel('Quadrature component');
    title('Time-varying phase of the complex baseband signal');

    % Limit the number of animation frames.
    numberOfAnimationFrames = min(500, length(time));
    animationIndex = round(linspace(1, length(time), ...
        numberOfAnimationFrames));

    for k = 1:length(animationIndex)

        index = animationIndex(k);

        currentI = real(basebandSignal(index));
        currentQ = imag(basebandSignal(index));

        addpoints(trajectory, currentI, currentQ);

        movingPoint.XData = currentI;
        movingPoint.YData = currentQ;

        phaseVector.UData = currentI;
        phaseVector.VData = currentQ;

        title(sprintf([ ...
            'Time = %.4f s, Phase = %.3f rad, ' ...
            'f_{inst} = %.3f Hz'], ...
            time(index), ...
            phaseExact(index), ...
            instantaneousFrequencyExact(index)));

        drawnow;
    end
end