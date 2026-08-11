function app = simple_human_vital_modeling(mode)
%SIMPLE_HUMAN_VITAL_MODELING Interactive far-field validation App.
%
%   SIMPLE_HUMAN_VITAL_MODELING opens a programmatic MATLAB App that
%   compares exact Euclidean and Chen far-field ranges for one point
%   scatterer on a translating rigid chest component. Respiration and
%   heartbeat parameters and the simulation duration can be changed from
%   the control panel.
%
%   APP = SIMPLE_HUMAN_VITAL_MODELING returns the App handles and the most
%   recently calculated results.
%
%   OUT = SIMPLE_HUMAN_VITAL_MODELING("selftest") runs the default model
%   and all numerical consistency checks without creating a user interface.
%
% Coordinate-frame hierarchy
% --------------------------
%   A: radar coordinates, with the monostatic radar at the origin
%   B: reference coordinates, parallel to A and centered at [R0;0;0]
%   C: target-local coordinates, with the same origin as B
%   D: rigid chest-component coordinates, translating relative to C
%
% The single point scatterer is fixed at the origin of D. Its position is
%
%   A_p(t) = A_T_B * B_T_C * C_T_D(t) * D_p.
%
% The vital translation is fixed along local +X_C = [1;0;0]. The fixed
% B-to-C rotation makes that local direction 45 degrees from the radar LOS,
% intentionally exposing the transverse-range term neglected by the
% far-field approximation. The Verification tab also sweeps this angle
% from 0 to 90 degrees.

    if nargin < 1
        mode = "app";
    end

    defaultParameters = makeDefaultParameters();

    if strcmpi(string(mode), "selftest")
        results = simulateVitalModel(defaultParameters);
        printSelfTestSummary(results);
        app = struct('Parameters', defaultParameters, ...
                     'Results', results, ...
                     'Validation', results.validation);
        return;
    elseif ~strcmpi(string(mode), "app")
        error('Unknown mode "%s". Use "app" or "selftest".', string(mode));
    end

    %% Build the App shell
    app = struct();
    app.Figure = uifigure( ...
        'Name', 'Human Vital-Sign Far-Field Validation', ...
        'Position', [70, 50, 1540, 900], ...
        'Color', [0.97, 0.97, 0.97]);

    mainGrid = uigridlayout(app.Figure, [1, 2]);
    mainGrid.ColumnWidth = {340, '1x'};
    mainGrid.Padding = [10, 10, 10, 10];
    mainGrid.ColumnSpacing = 10;

    leftGrid = uigridlayout(mainGrid, [3, 1]);
    leftGrid.Layout.Row = 1;
    leftGrid.Layout.Column = 1;
    leftGrid.RowHeight = {410, 135, '1x'};
    leftGrid.Padding = [0, 0, 0, 0];
    leftGrid.RowSpacing = 8;

    %% User parameter controls
    parameterPanel = uipanel(leftGrid, ...
        'Title', 'User parameters', ...
        'FontWeight', 'bold');
    parameterPanel.Layout.Row = 1;

    parameterGrid = uigridlayout(parameterPanel, [12, 2]);
    parameterGrid.ColumnWidth = {190, '1x'};
    parameterGrid.RowHeight = repmat({'fit'}, 1, 12);
    parameterGrid.Padding = [8, 8, 8, 8];

    app.InitialRangeField = addNumericControl(parameterGrid, 1, ...
        'Initial reference range R0 [m]', ...
        defaultParameters.initialRange, [0.1, 100]);

    app.SimulationTimeField = addNumericControl(parameterGrid, 2, ...
        'Simulation time [s]', ...
        defaultParameters.simulationTime, [1, 120]);

    addSectionLabel(parameterGrid, 3, 'Respiration');
    app.RespirationAmplitudeField = addNumericControl(parameterGrid, 4, ...
        'Peak displacement [mm]', ...
        1e3*defaultParameters.respirationAmplitude, [0, 20]);
    app.RespirationFrequencyField = addNumericControl(parameterGrid, 5, ...
        'Frequency [Hz]', ...
        defaultParameters.respirationFrequency, [0.05, 1]);
    app.RespirationPhaseField = addNumericControl(parameterGrid, 6, ...
        'Initial phase [deg]', ...
        rad2deg(defaultParameters.respirationPhase), [-360, 360]);

    addSectionLabel(parameterGrid, 7, 'Heartbeat');
    app.HeartbeatAmplitudeField = addNumericControl(parameterGrid, 8, ...
        'Peak displacement [mm]', ...
        1e3*defaultParameters.heartbeatAmplitude, [0, 5]);
    app.HeartbeatFrequencyField = addNumericControl(parameterGrid, 9, ...
        'Frequency [Hz]', ...
        defaultParameters.heartbeatFrequency, [0.5, 3]);
    app.HeartbeatPhaseField = addNumericControl(parameterGrid, 10, ...
        'Initial phase [deg]', ...
        rad2deg(defaultParameters.heartbeatPhase), [-360, 360]);

    app.UpdateButton = uibutton(parameterGrid, 'push', ...
        'Text', 'Update validation', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.22, 0.55, 0.82]);
    app.UpdateButton.Layout.Row = 11;
    app.UpdateButton.Layout.Column = 1;

    app.ResetButton = uibutton(parameterGrid, 'push', ...
        'Text', 'Reset defaults');
    app.ResetButton.Layout.Row = 11;
    app.ResetButton.Layout.Column = 2;

    app.AnimationButton = uibutton(parameterGrid, 'push', ...
        'Text', 'Animate scatterer', ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0.20, 0.65, 0.42]);
    app.AnimationButton.Layout.Row = 12;
    app.AnimationButton.Layout.Column = 1;

    app.StopAnimationButton = uibutton(parameterGrid, 'push', ...
        'Text', 'Stop animation', ...
        'Enable', 'off');
    app.StopAnimationButton.Layout.Row = 12;
    app.StopAnimationButton.Layout.Column = 2;

    %% Fixed modeling assumptions and run status
    settingsPanel = uipanel(leftGrid, ...
        'Title', 'Fixed model settings', ...
        'FontWeight', 'bold');
    settingsPanel.Layout.Row = 2;

    settingsGrid = uigridlayout(settingsPanel, [2, 2]);
    settingsGrid.ColumnWidth = {24, '1x'};
    settingsGrid.RowHeight = {24, '1x'};
    settingsGrid.Padding = [8, 6, 8, 6];

    app.StatusLamp = uilamp(settingsGrid, 'Color', [0.95, 0.65, 0.15]);
    app.StatusLamp.Layout.Row = 1;
    app.StatusLamp.Layout.Column = 1;

    app.StatusLabel = uilabel(settingsGrid, ...
        'Text', 'Waiting for calculation', ...
        'FontWeight', 'bold');
    app.StatusLabel.Layout.Row = 1;
    app.StatusLabel.Layout.Column = 2;

    fixedSettingsText = uitextarea(settingsGrid, ...
        'Editable', 'off', ...
        'Value', { ...
            sprintf('fc = %.1f GHz, monostatic, phase = +4*pi*DeltaR/lambda', ...
                    defaultParameters.centerFrequency/1e9), ...
            sprintf('u_C = [1 0 0]^T, vital/LOS angle = %.1f deg, fs = %.0f Hz', ...
                    defaultParameters.vitalToLosAngleDeg, ...
                    defaultParameters.samplingRate), ...
            'Unit reflectivity; propagation-amplitude loss is excluded'});
    fixedSettingsText.Layout.Row = 2;
    fixedSettingsText.Layout.Column = [1, 2];

    %% Numerical result monitor
    resultPanel = uipanel(leftGrid, ...
        'Title', 'Validation monitor', ...
        'FontWeight', 'bold');
    resultPanel.Layout.Row = 3;

    resultGrid = uigridlayout(resultPanel, [1, 1]);
    resultGrid.Padding = [5, 5, 5, 5];
    app.MetricTable = uitable(resultGrid, ...
        'ColumnName', {'Metric', 'Value'}, ...
        'ColumnWidth', {195, 110}, ...
        'RowName', {}, ...
        'Data', cell(0, 2));

    %% Plot tabs
    app.TabGroup = uitabgroup(mainGrid);
    app.TabGroup.Layout.Row = 1;
    app.TabGroup.Layout.Column = 2;

    rangeTab = uitab(app.TabGroup, 'Title', 'Range');
    rangeGrid = uigridlayout(rangeTab, [2, 2]);
    rangeGrid.Padding = [8, 8, 8, 8];
    app.MotionAxes = uiaxes(rangeGrid);
    app.RangeAxes = uiaxes(rangeGrid);
    app.RangeErrorAxes = uiaxes(rangeGrid);
    app.CalibratedErrorAxes = uiaxes(rangeGrid);

    phaseTab = uitab(app.TabGroup, 'Title', 'Phase / Baseband');
    phaseGrid = uigridlayout(phaseTab, [2, 2]);
    phaseGrid.Padding = [8, 8, 8, 8];
    app.PhaseAxes = uiaxes(phaseGrid);
    app.PhaseErrorAxes = uiaxes(phaseGrid);
    app.ComplexAxes = uiaxes(phaseGrid);
    app.SignalErrorAxes = uiaxes(phaseGrid);

    dopplerTab = uitab(app.TabGroup, ...
        'Title', 'Instantaneous frequency / Spectrogram');
    dopplerGrid = uigridlayout(dopplerTab, [2, 2]);
    dopplerGrid.Padding = [8, 8, 8, 8];
    app.InstantaneousFrequencyAxes = uiaxes(dopplerGrid);
    app.InstantaneousFrequencyErrorAxes = uiaxes(dopplerGrid);
    app.ExactSpectrogramAxes = uiaxes(dopplerGrid);
    app.FarFieldSpectrogramAxes = uiaxes(dopplerGrid);

    verificationTab = uitab(app.TabGroup, 'Title', 'Far-field verification');
    verificationGrid = uigridlayout(verificationTab, [2, 2]);
    verificationGrid.Padding = [8, 8, 8, 8];
    app.GeometryAxes = uiaxes(verificationGrid);
    app.DistanceSweepAxes = uiaxes(verificationGrid);
    app.AngleSweepAxes = uiaxes(verificationGrid);
    app.AsymptoticAxes = uiaxes(verificationGrid);

    animationTab = uitab(app.TabGroup, 'Title', 'Animation');
    animationGrid = uigridlayout(animationTab, [2, 2]);
    animationGrid.RowHeight = {'2x', '1x'};
    animationGrid.Padding = [8, 8, 8, 8];
    app.AnimationWholeAxes = uiaxes(animationGrid);
    app.AnimationCloseAxes = uiaxes(animationGrid);
    app.AnimationTimeAxes = uiaxes(animationGrid);
    app.AnimationTimeAxes.Layout.Row = 2;
    app.AnimationTimeAxes.Layout.Column = [1, 2];

    app.AnimationStopRequested = false;

    app.UpdateButton.ButtonPushedFcn = @updateApp;
    app.ResetButton.ButtonPushedFcn = @resetApp;
    app.AnimationButton.ButtonPushedFcn = @runScattererAnimation;
    app.StopAnimationButton.ButtonPushedFcn = @stopScattererAnimation;

    updateApp();

    %% Nested App callbacks
    function updateApp(~, ~)
        app.AnimationStopRequested = true;
        app.StatusLamp.Color = [0.95, 0.65, 0.15];
        app.StatusLabel.Text = 'Calculating...';
        app.UpdateButton.Enable = 'off';
        drawnow;

        try
            parameters = readParametersFromControls();
            results = simulateVitalModel(parameters);
            updateMetricTable(results);
            updatePlots(results);

            app.Parameters = parameters;
            app.Results = results;
            app.Validation = results.validation;
            app.AnimationStopRequested = false;
            app.StatusLamp.Color = [0.20, 0.70, 0.30];
            app.StatusLabel.Text = 'All numerical checks passed';
        catch exception
            app.StatusLamp.Color = [0.85, 0.20, 0.20];
            app.StatusLabel.Text = 'Validation failed';
            uialert(app.Figure, exception.message, ...
                'Vital-model validation error');
        end

        app.UpdateButton.Enable = 'on';
        drawnow limitrate;
    end

    function resetApp(~, ~)
        app.InitialRangeField.Value = defaultParameters.initialRange;
        app.SimulationTimeField.Value = defaultParameters.simulationTime;
        app.RespirationAmplitudeField.Value = ...
            1e3*defaultParameters.respirationAmplitude;
        app.RespirationFrequencyField.Value = ...
            defaultParameters.respirationFrequency;
        app.RespirationPhaseField.Value = ...
            rad2deg(defaultParameters.respirationPhase);
        app.HeartbeatAmplitudeField.Value = ...
            1e3*defaultParameters.heartbeatAmplitude;
        app.HeartbeatFrequencyField.Value = ...
            defaultParameters.heartbeatFrequency;
        app.HeartbeatPhaseField.Value = ...
            rad2deg(defaultParameters.heartbeatPhase);
        updateApp();
    end

    function parameters = readParametersFromControls()
        parameters = defaultParameters;
        parameters.initialRange = app.InitialRangeField.Value;
        parameters.simulationTime = app.SimulationTimeField.Value;
        parameters.respirationAmplitude = ...
            1e-3*app.RespirationAmplitudeField.Value;
        parameters.respirationFrequency = ...
            app.RespirationFrequencyField.Value;
        parameters.respirationPhase = ...
            deg2rad(app.RespirationPhaseField.Value);
        parameters.heartbeatAmplitude = ...
            1e-3*app.HeartbeatAmplitudeField.Value;
        parameters.heartbeatFrequency = ...
            app.HeartbeatFrequencyField.Value;
        parameters.heartbeatPhase = ...
            deg2rad(app.HeartbeatPhaseField.Value);
    end

    function updateMetricTable(results)
        metricData = { ...
            'Max |range error|', ...
                sprintf('%.6g um', 1e6*results.metrics.maximumRangeError); ...
            'Simulation time', ...
                sprintf('%.6g s', results.parameters.simulationTime); ...
            'RMS range error', ...
                sprintf('%.6g um', 1e6*results.metrics.rmsRangeError); ...
            'Max calibrated error', ...
                sprintf('%.6g um', 1e6*results.metrics.maximumCalibratedRangeError); ...
            'Max |phase error|', ...
                sprintf('%.6g mrad', 1e3*results.metrics.maximumPhaseError); ...
            'Max |phase error|', ...
                sprintf('%.6g deg', rad2deg(results.metrics.maximumPhaseError)); ...
            'RMS phase error', ...
                sprintf('%.6g mrad', 1e3*results.metrics.rmsPhaseError); ...
            'Max |instantaneous freq.|', ...
                sprintf('%.6g Hz', ...
                    results.metrics.maximumInstantaneousFrequency); ...
            'Max frequency error', ...
                sprintf('%.6g mHz', ...
                    1e3*results.metrics.maximumInstantaneousFrequencyError); ...
            'Max |baseband error|', ...
                sprintf('%.6g', results.metrics.maximumBasebandError); ...
            'RMS EVM', ...
                sprintf('%.6g %%', 100*results.metrics.rmsEvm); ...
            'Signal correlation', ...
                sprintf('%.12f', results.metrics.signalCorrelation); ...
            'Max displacement / R0', ...
                sprintf('%.6e', results.metrics.farFieldRatio); ...
            'Closed-form range check', ...
                sprintf('%.3e m', results.validation.maximumClosedFormError); ...
            'Transform-chain check', ...
                sprintf('%.3e m', results.validation.maximumTransformError); ...
            'Local +X direction check', ...
                sprintf('%.3e m', ...
                    results.validation.maximumLocalDirectionError)};
        app.MetricTable.Data = metricData;
    end

    function updatePlots(results)
        time = results.time;

        cla(app.MotionAxes);
        plot(app.MotionAxes, time, 1e3*results.motion.respiration, ...
            'LineWidth', 1.1, 'DisplayName', 'Respiration');
        hold(app.MotionAxes, 'on');
        plot(app.MotionAxes, time, 1e3*results.motion.heartbeat, ...
            'LineWidth', 1.1, 'DisplayName', 'Heartbeat');
        plot(app.MotionAxes, time, 1e3*results.motion.total, ...
            'k', 'LineWidth', 1.4, 'DisplayName', 'Total');
        hold(app.MotionAxes, 'off');
        styleTimeAxes(app.MotionAxes, 'Chest-component translation', ...
            'Displacement [mm]');
        legend(app.MotionAxes, 'Location', 'best');

        cla(app.RangeAxes);
        plot(app.RangeAxes, time, 1e3*results.range.exactVariation, ...
            'LineWidth', 1.4, 'DisplayName', 'Exact Euclidean');
        hold(app.RangeAxes, 'on');
        plot(app.RangeAxes, time, 1e3*results.range.farFieldVariation, ...
            '--', 'LineWidth', 1.3, 'DisplayName', 'Chen far-field');
        hold(app.RangeAxes, 'off');
        styleTimeAxes(app.RangeAxes, ...
            'Range variation relative to equilibrium R0', ...
            'Range variation [mm]');
        legend(app.RangeAxes, 'Location', 'best');

        cla(app.RangeErrorAxes);
        plot(app.RangeErrorAxes, time, 1e6*results.range.error, ...
            'Color', [0.78, 0.20, 0.18], 'LineWidth', 1.35);
        yline(app.RangeErrorAxes, 0, ':k');
        styleTimeAxes(app.RangeErrorAxes, ...
            'Raw geometric error: far-field - exact', ...
            'Range error [um]');

        cla(app.CalibratedErrorAxes);
        plot(app.CalibratedErrorAxes, time, ...
            1e6*results.range.calibratedError, ...
            'Color', [0.35, 0.20, 0.70], 'LineWidth', 1.35);
        yline(app.CalibratedErrorAxes, 0, ':k');
        styleTimeAxes(app.CalibratedErrorAxes, ...
            'Error after independent t = 0 range calibration', ...
            'Calibrated error [um]');

        cla(app.PhaseAxes);
        plot(app.PhaseAxes, time, results.phase.exact, ...
            'LineWidth', 1.35, 'DisplayName', 'Exact');
        hold(app.PhaseAxes, 'on');
        plot(app.PhaseAxes, time, results.phase.farField, ...
            '--', 'LineWidth', 1.25, 'DisplayName', 'Far-field');
        hold(app.PhaseAxes, 'off');
        styleTimeAxes(app.PhaseAxes, ...
            'Equilibrium-range-compensated phase', 'Phase [rad]');
        legend(app.PhaseAxes, 'Location', 'best');

        cla(app.PhaseErrorAxes);
        plot(app.PhaseErrorAxes, time, 1e3*results.phase.error, ...
            'Color', [0.78, 0.20, 0.18], 'LineWidth', 1.35);
        yline(app.PhaseErrorAxes, 0, ':k');
        styleTimeAxes(app.PhaseErrorAxes, ...
            'Phase error: far-field - exact', 'Phase error [mrad]');

        cla(app.ComplexAxes);
        plot(app.ComplexAxes, real(results.signal.exact), ...
            imag(results.signal.exact), 'LineWidth', 1.1, ...
            'DisplayName', 'Exact');
        hold(app.ComplexAxes, 'on');
        plot(app.ComplexAxes, real(results.signal.farField), ...
            imag(results.signal.farField), '--', 'LineWidth', 1.1, ...
            'DisplayName', 'Far-field');
        scatter(app.ComplexAxes, real(results.signal.exact(1)), ...
            imag(results.signal.exact(1)), 42, 'filled', ...
            'DisplayName', 't = 0');
        hold(app.ComplexAxes, 'off');
        axis(app.ComplexAxes, 'equal');
        grid(app.ComplexAxes, 'on');
        xlim(app.ComplexAxes, [-1.1, 1.1]);
        ylim(app.ComplexAxes, [-1.1, 1.1]);
        xlabel(app.ComplexAxes, 'In-phase');
        ylabel(app.ComplexAxes, 'Quadrature');
        title(app.ComplexAxes, 'Normalized complex baseband trajectory');
        legend(app.ComplexAxes, 'Location', 'best');

        cla(app.SignalErrorAxes);
        plot(app.SignalErrorAxes, time, abs(results.signal.error), ...
            'Color', [0.10, 0.50, 0.35], 'LineWidth', 1.35);
        styleTimeAxes(app.SignalErrorAxes, ...
            'Complex baseband error magnitude', '|s_{ff} - s_{exact}|');

        cla(app.InstantaneousFrequencyAxes);
        plot(app.InstantaneousFrequencyAxes, time, ...
            results.instantaneousFrequency.exact, ...
            'LineWidth', 1.35, 'DisplayName', 'Exact');
        hold(app.InstantaneousFrequencyAxes, 'on');
        plot(app.InstantaneousFrequencyAxes, time, ...
            results.instantaneousFrequency.farField, '--', ...
            'LineWidth', 1.25, 'DisplayName', 'Far-field');
        hold(app.InstantaneousFrequencyAxes, 'off');
        styleTimeAxes(app.InstantaneousFrequencyAxes, ...
            'Baseband instantaneous frequency', 'Frequency [Hz]');
        legend(app.InstantaneousFrequencyAxes, 'Location', 'best');

        cla(app.InstantaneousFrequencyErrorAxes);
        plot(app.InstantaneousFrequencyErrorAxes, time, ...
            1e3*results.instantaneousFrequency.error, ...
            'Color', [0.78, 0.20, 0.18], 'LineWidth', 1.35);
        yline(app.InstantaneousFrequencyErrorAxes, 0, ':k');
        styleTimeAxes(app.InstantaneousFrequencyErrorAxes, ...
            'Instantaneous-frequency error: far-field - exact', ...
            'Frequency error [mHz]');

        plotBasebandSpectrogram(app.ExactSpectrogramAxes, ...
            results.spectrogram.time, results.spectrogram.frequency, ...
            results.spectrogram.exactDb, results.spectrogram.displayLimit, ...
            'Exact Euclidean-range baseband spectrogram', time, ...
            results.instantaneousFrequency.exact);
        plotBasebandSpectrogram(app.FarFieldSpectrogramAxes, ...
            results.spectrogram.time, results.spectrogram.frequency, ...
            results.spectrogram.farFieldDb, ...
            results.spectrogram.displayLimit, ...
            'Chen far-field baseband spectrogram', time, ...
            results.instantaneousFrequency.farField);

        updateGeometryPlot(results);

        cla(app.DistanceSweepAxes);
        semilogx(app.DistanceSweepAxes, ...
            results.sweeps.distance.range, ...
            1e6*results.sweeps.distance.maximumError, ...
            'LineWidth', 1.4, 'DisplayName', 'Exact error');
        hold(app.DistanceSweepAxes, 'on');
        semilogx(app.DistanceSweepAxes, ...
            results.sweeps.distance.range, ...
            1e6*results.sweeps.distance.maximumSecondOrderError, ...
            '--', 'LineWidth', 1.2, 'DisplayName', 'Second order');
        xline(app.DistanceSweepAxes, results.parameters.initialRange, ':k', ...
            'Current R0', 'LabelVerticalAlignment', 'bottom');
        hold(app.DistanceSweepAxes, 'off');
        grid(app.DistanceSweepAxes, 'on');
        xlabel(app.DistanceSweepAxes, 'Reference range R0 [m]');
        ylabel(app.DistanceSweepAxes, 'Maximum |range error| [um]');
        title(app.DistanceSweepAxes, 'Far-field convergence with distance');
        legend(app.DistanceSweepAxes, 'Location', 'best');

        cla(app.AngleSweepAxes);
        plot(app.AngleSweepAxes, ...
            results.sweeps.angle.angleDeg, ...
            1e6*results.sweeps.angle.maximumError, ...
            'LineWidth', 1.4, 'DisplayName', 'Exact error');
        hold(app.AngleSweepAxes, 'on');
        plot(app.AngleSweepAxes, ...
            results.sweeps.angle.angleDeg, ...
            1e6*results.sweeps.angle.maximumSecondOrderError, ...
            '--', 'LineWidth', 1.2, 'DisplayName', 'Second order');
        xline(app.AngleSweepAxes, ...
            results.parameters.vitalToLosAngleDeg, ':k', 'Current angle', ...
            'LabelVerticalAlignment', 'bottom');
        hold(app.AngleSweepAxes, 'off');
        grid(app.AngleSweepAxes, 'on');
        xlabel(app.AngleSweepAxes, 'Vital-motion angle from LOS [deg]');
        ylabel(app.AngleSweepAxes, 'Maximum |range error| [um]');
        title(app.AngleSweepAxes, 'Neglected transverse component');
        legend(app.AngleSweepAxes, 'Location', 'best');

        cla(app.AsymptoticAxes);
        plot(app.AsymptoticAxes, time, 1e6*results.range.error, ...
            'LineWidth', 1.35, 'DisplayName', 'Exact error');
        hold(app.AsymptoticAxes, 'on');
        plot(app.AsymptoticAxes, time, ...
            1e6*results.range.secondOrderError, '--', ...
            'LineWidth', 1.2, 'DisplayName', ...
            '-|r_{perp}|^2/(2R0)');
        hold(app.AsymptoticAxes, 'off');
        styleTimeAxes(app.AsymptoticAxes, ...
            'Exact error versus far-field asymptote', 'Range error [um]');
        legend(app.AsymptoticAxes, 'Location', 'best');

        initializeAnimationPlots(results);
    end

    function runScattererAnimation(~, ~)
        if ~isfield(app, 'Results') || isempty(app.Results)
            return;
        end

        results = app.Results;
        app.AnimationStopRequested = false;
        app.AnimationButton.Enable = 'off';
        app.StopAnimationButton.Enable = 'on';
        app.UpdateButton.Enable = 'off';
        app.ResetButton.Enable = 'off';
        app.TabGroup.SelectedTab = animationTab;

        maximumFrames = 240;
        numberOfSamples = numel(results.time);
        numberOfFrames = min(maximumFrames, numberOfSamples);
        frameIndices = unique(round(linspace(1, numberOfSamples, ...
            numberOfFrames)));
        playbackDuration = min(10, max(5, 0.35*results.parameters.simulationTime));
        framePause = playbackDuration/max(1, numel(frameIndices)-1);

        initializeAnimationPlots(results);
        app.StatusLamp.Color = [0.20, 0.55, 0.90];

        for frameIndex = frameIndices
            if app.AnimationStopRequested || ~isvalid(app.Figure)
                break;
            end

            updateAnimationFrame(results, frameIndex);
            app.StatusLabel.Text = sprintf('Animation: t = %.3f / %.3f s', ...
                results.time(frameIndex), results.time(end));
            drawnow;
            pause(framePause);
        end

        if isvalid(app.Figure)
            app.AnimationButton.Enable = 'on';
            app.StopAnimationButton.Enable = 'off';
            app.UpdateButton.Enable = 'on';
            app.ResetButton.Enable = 'on';
            app.AnimationStopRequested = false;
            app.StatusLamp.Color = [0.20, 0.70, 0.30];
            app.StatusLabel.Text = 'All numerical checks passed';
        end
    end

    function stopScattererAnimation(~, ~)
        app.AnimationStopRequested = true;
        app.StatusLabel.Text = 'Stopping animation...';
    end

    function initializeAnimationPlots(results)
        wholeAxes = app.AnimationWholeAxes;
        closeAxes = app.AnimationCloseAxes;
        timeAxes = app.AnimationTimeAxes;

        radarPosition = results.geometry.radarPosition_A;
        referenceOrigin = results.geometry.referenceOrigin_A;
        pointPath_A = results.geometry.pointPosition_A;

        cla(wholeAxes);
        plot3(wholeAxes, ...
            [radarPosition(1), referenceOrigin(1)], ...
            [radarPosition(2), referenceOrigin(2)], ...
            [radarPosition(3), referenceOrigin(3)], ...
            'k--', 'LineWidth', 1.1, 'DisplayName', 'Reference LOS');
        hold(wholeAxes, 'on');
        plot3(wholeAxes, pointPath_A(1,:), pointPath_A(2,:), ...
            pointPath_A(3,:), 'Color', [0.75, 0.75, 0.75], ...
            'LineWidth', 1.1, 'DisplayName', 'Full scatterer path');
        scatter3(wholeAxes, radarPosition(1), radarPosition(2), ...
            radarPosition(3), 70, 'k', 'filled', 'DisplayName', 'Radar');
        scatter3(wholeAxes, referenceOrigin(1), referenceOrigin(2), ...
            referenceOrigin(3), 65, [0.1, 0.4, 0.8], 'filled', ...
            'DisplayName', 'B/C origin');
        app.AnimationWholePoint = scatter3(wholeAxes, ...
            pointPath_A(1,1), pointPath_A(2,1), pointPath_A(3,1), ...
            85, [0.85, 0.20, 0.15], 'filled', ...
            'DisplayName', 'Moving point');
        app.AnimationExactLine = plot3(wholeAxes, ...
            [radarPosition(1), pointPath_A(1,1)], ...
            [radarPosition(2), pointPath_A(2,1)], ...
            [radarPosition(3), pointPath_A(3,1)], ...
            'Color', [0.85, 0.20, 0.15], 'LineWidth', 1.4, ...
            'DisplayName', 'Exact range');
        hold(wholeAxes, 'off');
        grid(wholeAxes, 'on');
        axis(wholeAxes, 'equal');
        view(wholeAxes, 30, 22);
        xlabel(wholeAxes, 'X_A [m]');
        ylabel(wholeAxes, 'Y_A [m]');
        zlabel(wholeAxes, 'Z_A [m]');
        title(wholeAxes, 'Radar coordinates: complete geometry');
        legend(wholeAxes, 'Location', 'best');

        pointPath_C_mm = 1e3*results.geometry.pointPosition_C;
        closeLimit = max(1, 1.20*max(abs(pointPath_C_mm(1,:))));
        transverseLimit = max(0.25, 0.12*closeLimit);

        cla(closeAxes);
        plot3(closeAxes, pointPath_C_mm(1,:), pointPath_C_mm(2,:), ...
            pointPath_C_mm(3,:), 'Color', [0.70, 0.70, 0.70], ...
            'LineWidth', 1.3, 'DisplayName', 'Scatterer path');
        hold(closeAxes, 'on');
        quiver3(closeAxes, 0, 0, 0, 0.75*closeLimit, 0, 0, 0, ...
            'Color', [0.10, 0.55, 0.30], 'LineWidth', 1.5, ...
            'DisplayName', '+X_C vital direction');
        scatter3(closeAxes, 0, 0, 0, 55, [0.1, 0.4, 0.8], 'filled', ...
            'DisplayName', 'C equilibrium origin');
        app.AnimationClosePoint = scatter3(closeAxes, ...
            pointPath_C_mm(1,1), pointPath_C_mm(2,1), ...
            pointPath_C_mm(3,1), 95, [0.85, 0.20, 0.15], 'filled', ...
            'DisplayName', 'Point scatterer');
        hold(closeAxes, 'off');
        grid(closeAxes, 'on');
        axis(closeAxes, 'equal');
        view(closeAxes, 28, 20);
        xlim(closeAxes, [-closeLimit, closeLimit]);
        ylim(closeAxes, [-transverseLimit, transverseLimit]);
        zlim(closeAxes, [-transverseLimit, transverseLimit]);
        xlabel(closeAxes, 'X_C [mm]');
        ylabel(closeAxes, 'Y_C [mm]');
        zlabel(closeAxes, 'Z_C [mm]');
        title(closeAxes, 'Local coordinates: magnified chest motion');
        legend(closeAxes, 'Location', 'best');

        cla(timeAxes);
        plot(timeAxes, results.time, 1e3*results.motion.total, ...
            'k', 'LineWidth', 1.25, 'DisplayName', 'Total vital motion');
        hold(timeAxes, 'on');
        app.AnimationTimeMarker = scatter(timeAxes, results.time(1), ...
            1e3*results.motion.total(1), 70, [0.85, 0.20, 0.15], ...
            'filled', 'DisplayName', 'Current point');
        app.AnimationTimeLine = xline(timeAxes, results.time(1), '--', ...
            'Color', [0.20, 0.45, 0.80], 'LineWidth', 1.1, ...
            'DisplayName', 'Current time');
        hold(timeAxes, 'off');
        styleTimeAxes(timeAxes, 'Animation timeline', ...
            'Local +X_C displacement [mm]');
        legend(timeAxes, 'Location', 'best');

        updateAnimationFrame(results, 1);
    end

    function updateAnimationFrame(results, timeIndex)
        point_A = results.geometry.pointPosition_A(:,timeIndex);
        point_C_mm = 1e3*results.geometry.pointPosition_C(:,timeIndex);
        radarPosition = results.geometry.radarPosition_A;

        app.AnimationWholePoint.XData = point_A(1);
        app.AnimationWholePoint.YData = point_A(2);
        app.AnimationWholePoint.ZData = point_A(3);
        app.AnimationExactLine.XData = [radarPosition(1), point_A(1)];
        app.AnimationExactLine.YData = [radarPosition(2), point_A(2)];
        app.AnimationExactLine.ZData = [radarPosition(3), point_A(3)];

        app.AnimationClosePoint.XData = point_C_mm(1);
        app.AnimationClosePoint.YData = point_C_mm(2);
        app.AnimationClosePoint.ZData = point_C_mm(3);

        app.AnimationTimeMarker.XData = results.time(timeIndex);
        app.AnimationTimeMarker.YData = 1e3*results.motion.total(timeIndex);
        app.AnimationTimeLine.Value = results.time(timeIndex);
        title(app.AnimationCloseAxes, sprintf( ...
            'Local coordinates: p_C = [%.4f, %.4f, %.4f]^T mm', ...
            point_C_mm(1), point_C_mm(2), point_C_mm(3)));
    end

    function updateGeometryPlot(results)
        axesHandle = app.GeometryAxes;
        cla(axesHandle);

        [~, pointIndex] = max(abs(results.motion.total));
        radarPosition = results.geometry.radarPosition_A;
        referenceOrigin = results.geometry.referenceOrigin_A;
        pointPosition = results.geometry.pointPosition_A(:, pointIndex);
        motionDirection = results.geometry.vitalDirection_A;
        directionLength = max(0.15*results.parameters.initialRange, 0.05);

        plot3(axesHandle, ...
            [radarPosition(1), referenceOrigin(1)], ...
            [radarPosition(2), referenceOrigin(2)], ...
            [radarPosition(3), referenceOrigin(3)], ...
            'k--', 'LineWidth', 1.15, 'DisplayName', 'Reference LOS');
        hold(axesHandle, 'on');
        scatter3(axesHandle, radarPosition(1), radarPosition(2), ...
            radarPosition(3), 75, 'k', 'filled', 'DisplayName', 'Radar');
        scatter3(axesHandle, referenceOrigin(1), referenceOrigin(2), ...
            referenceOrigin(3), 70, [0.1, 0.4, 0.8], 'filled', ...
            'DisplayName', 'B/C origin');
        plot3(axesHandle, results.geometry.pointPosition_A(1,:), ...
            results.geometry.pointPosition_A(2,:), ...
            results.geometry.pointPosition_A(3,:), ...
            'Color', [0.82, 0.25, 0.18], 'LineWidth', 1.5, ...
            'DisplayName', 'Scatterer path');
        scatter3(axesHandle, pointPosition(1), pointPosition(2), ...
            pointPosition(3), 65, [0.82, 0.25, 0.18], 'filled', ...
            'DisplayName', 'Point scatterer');
        quiver3(axesHandle, referenceOrigin(1), referenceOrigin(2), ...
            referenceOrigin(3), directionLength*motionDirection(1), ...
            directionLength*motionDirection(2), ...
            directionLength*motionDirection(3), 0, ...
            'Color', [0.10, 0.55, 0.30], 'LineWidth', 1.6, ...
            'MaxHeadSize', 0.5, ...
            'DisplayName', 'Local +X_C vital direction');
        hold(axesHandle, 'off');

        grid(axesHandle, 'on');
        axis(axesHandle, 'equal');
        view(axesHandle, 30, 22);
        xlabel(axesHandle, 'X_A [m]');
        ylabel(axesHandle, 'Y_A [m]');
        zlabel(axesHandle, 'Z_A [m]');
        title(axesHandle, sprintf('A-B-C-D geometry (motion angle %.1f deg)', ...
            results.parameters.vitalToLosAngleDeg));
        legend(axesHandle, 'Location', 'best');
    end
end

%% UI helpers
function field = addNumericControl(parent, row, labelText, value, limits)
    label = uilabel(parent, 'Text', labelText);
    label.Layout.Row = row;
    label.Layout.Column = 1;

    field = uieditfield(parent, 'numeric', ...
        'Value', value, ...
        'Limits', limits, ...
        'LowerLimitInclusive', 'on', ...
        'UpperLimitInclusive', 'on');
    field.Layout.Row = row;
    field.Layout.Column = 2;
end

function addSectionLabel(parent, row, labelText)
    label = uilabel(parent, ...
        'Text', labelText, ...
        'FontWeight', 'bold', ...
        'FontColor', [0.15, 0.35, 0.60]);
    label.Layout.Row = row;
    label.Layout.Column = [1, 2];
end

function styleTimeAxes(axesHandle, titleText, yLabelText)
    grid(axesHandle, 'on');
    xlabel(axesHandle, 'Time [s]');
    ylabel(axesHandle, yLabelText);
    title(axesHandle, titleText);
end

function plotBasebandSpectrogram(axesHandle, time, frequency, ...
        magnitudeDb, frequencyLimit, titleText, signalTime, ...
        instantaneousFrequency)
    cla(axesHandle);
    imagesc(axesHandle, time, frequency, magnitudeDb);
    hold(axesHandle, 'on');
    plot(axesHandle, signalTime, instantaneousFrequency, 'w-', ...
        'LineWidth', 1.15);
    hold(axesHandle, 'off');
    axis(axesHandle, 'xy');
    ylim(axesHandle, [-frequencyLimit, frequencyLimit]);
    clim(axesHandle, [-60, 0]);
    colormap(axesHandle, turbo(256));
    colorbar(axesHandle);
    xlabel(axesHandle, 'Time [s]');
    ylabel(axesHandle, 'Frequency [Hz]');
    title(axesHandle, titleText);
end

%% Model and validation
function parameters = makeDefaultParameters()
    parameters.initialRange = 2.0;
    parameters.simulationTime = 20;
    parameters.respirationAmplitude = 5e-3;
    parameters.respirationFrequency = 0.25;
    parameters.respirationPhase = deg2rad(0);
    parameters.heartbeatAmplitude = 0.5e-3;
    parameters.heartbeatFrequency = 1.2;
    parameters.heartbeatPhase = deg2rad(0);

    parameters.speedOfLight = 299792458;
    parameters.centerFrequency = 24e9;
    parameters.samplingRate = 200;
    parameters.vitalToLosAngleDeg = 0;
end

function results = simulateVitalModel(parameters)
    validateModelParameters(parameters);

    wavelength = parameters.speedOfLight / parameters.centerFrequency;
    phasePerMeter = 4*pi/wavelength;

    numberOfTimeIntervals = max(1, ...
        round(parameters.simulationTime*parameters.samplingRate));
    time = linspace(0, parameters.simulationTime, numberOfTimeIntervals+1);
    numberOfTimes = numel(time);

    respiration = parameters.respirationAmplitude * ...
        sin(2*pi*parameters.respirationFrequency*time + ...
            parameters.respirationPhase);
    heartbeat = parameters.heartbeatAmplitude * ...
        sin(2*pi*parameters.heartbeatFrequency*time + ...
            parameters.heartbeatPhase);
    totalDisplacement = respiration + heartbeat;
    respirationVelocity = 2*pi*parameters.respirationFrequency * ...
        parameters.respirationAmplitude .* ...
        cos(2*pi*parameters.respirationFrequency*time + ...
            parameters.respirationPhase);
    heartbeatVelocity = 2*pi*parameters.heartbeatFrequency * ...
        parameters.heartbeatAmplitude .* ...
        cos(2*pi*parameters.heartbeatFrequency*time + ...
            parameters.heartbeatPhase);
    totalVelocity = respirationVelocity + heartbeatVelocity;

    % A: radar, B: reference, C: local, D: chest component.
    radarPosition_A = [0; 0; 0];
    referenceOrigin_A = [parameters.initialRange; 0; 0];
    transform_A_B = makeTransform(eye(3), referenceOrigin_A);

    rotation_B_C = rotationZ(deg2rad(parameters.vitalToLosAngleDeg));
    transform_B_C = makeTransform(rotation_B_C, zeros(3,1));

    vitalDirection_C = [1; 0; 0];
    vitalDirection_A = rotation_B_C*vitalDirection_C;
    pointPosition_D = [0; 0; 0];
    pointPositionHomogeneous_D = [pointPosition_D; 1];

    referenceVector_A = referenceOrigin_A - radarPosition_A;
    referenceRange = norm(referenceVector_A);
    lineOfSight_A = referenceVector_A/referenceRange;
    projectionFactor = dot(lineOfSight_A, vitalDirection_A);

    pointPosition_A = zeros(3, numberOfTimes);
    pointPosition_C = zeros(3, numberOfTimes);
    transform_A_D = zeros(4, 4, numberOfTimes);
    exactRange = zeros(1, numberOfTimes);
    farFieldRange = zeros(1, numberOfTimes);
    expandedPositionError = zeros(1, numberOfTimes);

    for timeIndex = 1:numberOfTimes
        componentOrigin_C = totalDisplacement(timeIndex)*vitalDirection_C;
        currentPointPosition_C = componentOrigin_C + pointPosition_D;
        transform_C_D = makeTransform(eye(3), componentOrigin_C);
        currentTransform_A_D = transform_A_B*transform_B_C*transform_C_D;
        currentPointHomogeneous_A = ...
            currentTransform_A_D*pointPositionHomogeneous_D;
        currentPointPosition_A = currentPointHomogeneous_A(1:3);

        pointPosition_A(:,timeIndex) = currentPointPosition_A;
        pointPosition_C(:,timeIndex) = currentPointPosition_C;
        transform_A_D(:,:,timeIndex) = currentTransform_A_D;
        exactRange(timeIndex) = norm(currentPointPosition_A-radarPosition_A);

        relativePointPosition_A = currentPointPosition_A-referenceOrigin_A;
        farFieldRange(timeIndex) = referenceRange + ...
            dot(lineOfSight_A, relativePointPosition_A);

        expandedPointPosition_A = referenceOrigin_A + ...
            rotation_B_C*(componentOrigin_C + pointPosition_D);
        expandedPositionError(timeIndex) = ...
            norm(expandedPointPosition_A-currentPointPosition_A);
    end

    closedFormFarFieldRange = referenceRange + ...
        projectionFactor*totalDisplacement;
    exactVariation = exactRange-referenceRange;
    farFieldVariation = farFieldRange-referenceRange;
    rangeError = farFieldRange-exactRange;
    calibratedRangeError = rangeError-rangeError(1);

    transverseFractionSquared = max(0, 1-projectionFactor^2);
    transverseRangeSquared = ...
        transverseFractionSquared*totalDisplacement.^2;
    stableRangeError = -transverseRangeSquared ./ ...
        (exactRange+farFieldRange);
    secondOrderRangeError = -transverseRangeSquared/(2*referenceRange);

    exactPhase = phasePerMeter*exactVariation;
    farFieldPhase = phasePerMeter*farFieldVariation;
    phaseError = farFieldPhase-exactPhase;

    exactSignal = exp(1j*exactPhase);
    farFieldSignal = exp(1j*farFieldPhase);
    signalError = farFieldSignal-exactSignal;

    pointVelocity_A = vitalDirection_A*totalVelocity;
    radarToPointVector_A = pointPosition_A-radarPosition_A;
    exactRangeRate = sum(radarToPointVector_A.*pointVelocity_A, 1) ./ ...
        exactRange;
    farFieldRangeRate = projectionFactor*totalVelocity;
    exactInstantaneousFrequency = (2/wavelength)*exactRangeRate;
    farFieldInstantaneousFrequency = (2/wavelength)*farFieldRangeRate;
    instantaneousFrequencyError = farFieldInstantaneousFrequency - ...
        exactInstantaneousFrequency;

    maximumVitalFrequency = max(parameters.respirationFrequency, ...
                                parameters.heartbeatFrequency);
    maximumDisplayedInstantaneousFrequency = max(abs([ ...
        exactInstantaneousFrequency, farFieldInstantaneousFrequency]));
    spectrogramDisplayLimit = min(0.95*parameters.samplingRate/2, ...
        max(5, 1.5*maximumDisplayedInstantaneousFrequency + ...
            2*maximumVitalFrequency));
    spectrogramResults = calculateBasebandSpectrogram( ...
        exactSignal, farFieldSignal, parameters.samplingRate, ...
        spectrogramDisplayLimit);

    metrics.maximumRangeError = max(abs(rangeError));
    metrics.rmsRangeError = sqrt(mean(abs(rangeError).^2));
    metrics.maximumCalibratedRangeError = ...
        max(abs(calibratedRangeError));
    metrics.maximumPhaseError = max(abs(phaseError));
    metrics.rmsPhaseError = sqrt(mean(abs(phaseError).^2));
    metrics.maximumInstantaneousFrequency = ...
        maximumDisplayedInstantaneousFrequency;
    metrics.maximumInstantaneousFrequencyError = ...
        max(abs(instantaneousFrequencyError));
    metrics.rmsInstantaneousFrequencyError = ...
        sqrt(mean(abs(instantaneousFrequencyError).^2));
    metrics.maximumBasebandError = max(abs(signalError));
    metrics.rmsEvm = sqrt(mean(abs(signalError).^2)) / ...
        sqrt(mean(abs(exactSignal).^2));
    metrics.signalCorrelation = abs(sum(conj(exactSignal).*farFieldSignal)) / ...
        (norm(exactSignal)*norm(farFieldSignal));
    metrics.farFieldRatio = max(abs(totalDisplacement))/referenceRange;

    distanceSweep = calculateDistanceSweep(parameters, totalDisplacement, ...
        projectionFactor);
    angleSweep = calculateAngleSweep(parameters, totalDisplacement);

    validation.maximumTransformError = max(expandedPositionError);
    recoveredPointPosition_C = rotation_B_C.' * ...
        (pointPosition_A-referenceOrigin_A);
    expectedPointPosition_C = vitalDirection_C*totalDisplacement;
    validation.maximumLocalDirectionError = max(vecnorm( ...
        recoveredPointPosition_C-expectedPointPosition_C, 2, 1));
    validation.maximumClosedFormError = ...
        max(abs(farFieldRange-closedFormFarFieldRange));
    validation.maximumStableRangeIdentityError = ...
        max(abs(rangeError-stableRangeError));
    validation.maximumPhaseRangeIdentityError = ...
        max(abs(phaseError-phasePerMeter*rangeError));
    validation.maximumBasebandIdentityError = max(abs( ...
        abs(signalError)-2*abs(sin(phaseError/2))));
    validation.maximumFrequencyRangeRateIdentityError = max(abs( ...
        instantaneousFrequencyError - ...
        (2/wavelength)*(farFieldRangeRate-exactRangeRate)));
    validation.spectrogramIsFinite = ...
        all(isfinite(spectrogramResults.exactDb(:))) && ...
        all(isfinite(spectrogramResults.farFieldDb(:)));
    validation.maximumParallelMotionError = ...
        calculateParallelMotionError(referenceRange, totalDisplacement);

    validation.tolerances.transform = 1e-12;
    validation.tolerances.localDirection = 1e-12;
    validation.tolerances.closedForm = 1e-12;
    validation.tolerances.rangeIdentity = 1e-12;
    validation.tolerances.phaseIdentity = 1e-12;
    validation.tolerances.basebandIdentity = 1e-12;
    validation.tolerances.frequencyIdentity = 1e-12;
    validation.tolerances.parallelMotion = 1e-12;
    validation.allPassed = ...
        validation.maximumTransformError < validation.tolerances.transform && ...
        validation.maximumLocalDirectionError < ...
            validation.tolerances.localDirection && ...
        validation.maximumClosedFormError < validation.tolerances.closedForm && ...
        validation.maximumStableRangeIdentityError < ...
            validation.tolerances.rangeIdentity && ...
        validation.maximumPhaseRangeIdentityError < ...
            validation.tolerances.phaseIdentity && ...
        validation.maximumBasebandIdentityError < ...
            validation.tolerances.basebandIdentity && ...
        validation.maximumFrequencyRangeRateIdentityError < ...
            validation.tolerances.frequencyIdentity && ...
        validation.spectrogramIsFinite && ...
        validation.maximumParallelMotionError < ...
            validation.tolerances.parallelMotion;

    assert(validation.allPassed, ...
        'At least one geometry, range, phase, or baseband check failed.');

    results.parameters = parameters;
    results.time = time;
    results.wavelength = wavelength;
    results.phasePerMeter = phasePerMeter;
    results.motion.respiration = respiration;
    results.motion.heartbeat = heartbeat;
    results.motion.total = totalDisplacement;
    results.motion.respirationVelocity = respirationVelocity;
    results.motion.heartbeatVelocity = heartbeatVelocity;
    results.motion.totalVelocity = totalVelocity;
    results.geometry.radarPosition_A = radarPosition_A;
    results.geometry.referenceOrigin_A = referenceOrigin_A;
    results.geometry.lineOfSight_A = lineOfSight_A;
    results.geometry.vitalDirection_C = vitalDirection_C;
    results.geometry.vitalDirection_A = vitalDirection_A;
    results.geometry.pointPosition_D = pointPosition_D;
    results.geometry.pointPosition_C = pointPosition_C;
    results.geometry.pointPosition_A = pointPosition_A;
    results.geometry.transform_A_B = transform_A_B;
    results.geometry.transform_B_C = transform_B_C;
    results.geometry.transform_A_D = transform_A_D;
    results.geometry.projectionFactor = projectionFactor;
    results.range.exact = exactRange;
    results.range.farField = farFieldRange;
    results.range.closedFormFarField = closedFormFarFieldRange;
    results.range.exactVariation = exactVariation;
    results.range.farFieldVariation = farFieldVariation;
    results.range.error = rangeError;
    results.range.calibratedError = calibratedRangeError;
    results.range.stableIdentityError = stableRangeError;
    results.range.secondOrderError = secondOrderRangeError;
    results.phase.exact = exactPhase;
    results.phase.farField = farFieldPhase;
    results.phase.error = phaseError;
    results.phase.wrappedError = angle(farFieldSignal.*conj(exactSignal));
    results.signal.exact = exactSignal;
    results.signal.farField = farFieldSignal;
    results.signal.error = signalError;
    results.instantaneousFrequency.exactRangeRate = exactRangeRate;
    results.instantaneousFrequency.farFieldRangeRate = farFieldRangeRate;
    results.instantaneousFrequency.exact = exactInstantaneousFrequency;
    results.instantaneousFrequency.farField = farFieldInstantaneousFrequency;
    results.instantaneousFrequency.error = instantaneousFrequencyError;
    results.spectrogram = spectrogramResults;
    results.metrics = metrics;
    results.sweeps.distance = distanceSweep;
    results.sweeps.angle = angleSweep;
    results.validation = validation;
end

function results = calculateBasebandSpectrogram(exactSignal, ...
        farFieldSignal, samplingRate, displayLimit)
    numberOfSamples = numel(exactSignal);
    windowLength = min(1024, max(32, floor(numberOfSamples/4)));
    overlapLength = min(windowLength-1, round(0.80*windowLength));
    fftLength = max(512, 2^nextpow2(2*windowLength));

    windowIndex = (0:windowLength-1).';
    analysisWindow = 0.54 - 0.46*cos( ...
        2*pi*windowIndex/max(1, windowLength-1));

    [exactStft, frequency, spectrogramTime] = spectrogram( ...
        exactSignal, analysisWindow, overlapLength, fftLength, ...
        samplingRate, 'centered');
    farFieldStft = spectrogram(farFieldSignal, analysisWindow, ...
        overlapLength, fftLength, samplingRate, 'centered');

    commonMagnitudeReference = max([abs(exactStft(:)); ...
                                    abs(farFieldStft(:))]);
    if commonMagnitudeReference == 0
        commonMagnitudeReference = 1;
    end

    results.time = spectrogramTime;
    results.frequency = frequency;
    results.exactStft = exactStft;
    results.farFieldStft = farFieldStft;
    results.exactDb = 20*log10(abs(exactStft)/commonMagnitudeReference + eps);
    results.farFieldDb = ...
        20*log10(abs(farFieldStft)/commonMagnitudeReference + eps);
    results.displayLimit = displayLimit;
    results.windowLength = windowLength;
    results.overlapLength = overlapLength;
    results.fftLength = fftLength;
end

function sweep = calculateDistanceSweep(parameters, displacement, projectionFactor)
    minimumRange = max(0.1, parameters.initialRange/10);
    maximumRange = min(1000, max(10*minimumRange, ...
        10*parameters.initialRange));
    rangeValues = logspace(log10(minimumRange), log10(maximumRange), 90);

    maximumError = zeros(size(rangeValues));
    maximumSecondOrderError = zeros(size(rangeValues));
    transverseFractionSquared = max(0, 1-projectionFactor^2);

    for rangeIndex = 1:numel(rangeValues)
        currentRange = rangeValues(rangeIndex);
        exactRange = sqrt(currentRange^2 + ...
            2*currentRange*projectionFactor*displacement + displacement.^2);
        farFieldRange = currentRange + projectionFactor*displacement;
        currentError = farFieldRange-exactRange;
        secondOrderError = -transverseFractionSquared*displacement.^2 / ...
            (2*currentRange);
        maximumError(rangeIndex) = max(abs(currentError));
        maximumSecondOrderError(rangeIndex) = max(abs(secondOrderError));
    end

    sweep.range = rangeValues;
    sweep.maximumError = maximumError;
    sweep.maximumSecondOrderError = maximumSecondOrderError;
end

function sweep = calculateAngleSweep(parameters, displacement)
    angleDeg = linspace(0, 90, 91);
    maximumError = zeros(size(angleDeg));
    maximumSecondOrderError = zeros(size(angleDeg));
    referenceRange = parameters.initialRange;

    for angleIndex = 1:numel(angleDeg)
        projectionFactor = cosd(angleDeg(angleIndex));
        exactRange = sqrt(referenceRange^2 + ...
            2*referenceRange*projectionFactor*displacement + displacement.^2);
        farFieldRange = referenceRange + projectionFactor*displacement;
        currentError = farFieldRange-exactRange;
        secondOrderError = -(1-projectionFactor^2)*displacement.^2 / ...
            (2*referenceRange);
        maximumError(angleIndex) = max(abs(currentError));
        maximumSecondOrderError(angleIndex) = max(abs(secondOrderError));
    end

    sweep.angleDeg = angleDeg;
    sweep.maximumError = maximumError;
    sweep.maximumSecondOrderError = maximumSecondOrderError;
end

function maximumError = calculateParallelMotionError(referenceRange, displacement)
    exactRange = sqrt(referenceRange^2 + ...
        2*referenceRange*displacement + displacement.^2);
    farFieldRange = referenceRange+displacement;
    maximumError = max(abs(farFieldRange-exactRange));
end

function validateModelParameters(parameters)
    validateattributes(parameters.initialRange, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(parameters.simulationTime, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(parameters.respirationAmplitude, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(parameters.heartbeatAmplitude, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(parameters.respirationFrequency, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});
    validateattributes(parameters.heartbeatFrequency, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'positive'});

    maximumDisplacement = parameters.respirationAmplitude + ...
                          parameters.heartbeatAmplitude;
    assert(maximumDisplacement < parameters.initialRange, ...
        'The total vital displacement must be smaller than R0.');
end

function transform = makeTransform(rotation, translation)
    transform = [rotation, translation(:); 0, 0, 0, 1];
end

function rotation = rotationZ(angle)
    rotation = [cos(angle), -sin(angle), 0; ...
                sin(angle),  cos(angle), 0; ...
                0,           0,          1];
end

function printSelfTestSummary(results)
    fprintf('============================================================\n');
    fprintf('Human vital-sign far-field model self-test\n');
    fprintf('============================================================\n');
    fprintf('Initial reference range       : %.6f m\n', ...
        results.parameters.initialRange);
    fprintf('Simulation time               : %.6f s\n', ...
        results.parameters.simulationTime);
    fprintf('Vital-motion/LOS angle        : %.3f deg\n', ...
        results.parameters.vitalToLosAngleDeg);
    fprintf('Maximum absolute range error  : %.9e m\n', ...
        results.metrics.maximumRangeError);
    fprintf('RMS range error               : %.9e m\n', ...
        results.metrics.rmsRangeError);
    fprintf('Maximum absolute phase error  : %.9e rad\n', ...
        results.metrics.maximumPhaseError);
    fprintf('Maximum instantaneous freq.   : %.9e Hz\n', ...
        results.metrics.maximumInstantaneousFrequency);
    fprintf('Maximum inst. freq. error     : %.9e Hz\n', ...
        results.metrics.maximumInstantaneousFrequencyError);
    fprintf('RMS baseband EVM              : %.9e\n', ...
        results.metrics.rmsEvm);
    fprintf('Closed-form range check       : %.3e m\n', ...
        results.validation.maximumClosedFormError);
    fprintf('Transform-chain check         : %.3e m\n', ...
        results.validation.maximumTransformError);
    fprintf('Local +X direction check      : %.3e m\n', ...
        results.validation.maximumLocalDirectionError);
    fprintf('Frequency/range-rate check    : %.3e Hz\n', ...
        results.validation.maximumFrequencyRangeRateIdentityError);
    fprintf('All checks passed             : %s\n', ...
        string(results.validation.allPassed));
end
