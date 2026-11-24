%% EXPERIMENT 01: ECG NOISE REMOVAL AND QRS DETECTION
% Course: EC4414 / 19BM414 | CO1
% Educational demonstration using a reproducible synthetic ECG-like signal.
% This is NOT a clinical-grade ECG detector.
clear; close all; clc; rng(11);

fs = 250;                 % Sampling frequency (Hz)
duration = 10;            % Signal duration (s)
t = (0:1/fs:duration-1/fs)';
cleanECG = zeros(size(t));
rPeakTimes = 0.7:0.85:9.2;

% Construct a simple P-QRS-T morphology around each R peak.
for k = 1:numel(rPeakTimes)
    r = rPeakTimes(k);
    cleanECG = cleanECG ...
        + 0.10*exp(-((t-(r-0.20))/0.045).^2) ... % P wave
        - 0.15*exp(-((t-(r-0.045))/0.018).^2) ... % Q wave
        + 1.00*exp(-((t-r)/0.020).^2) ...          % R wave
        - 0.25*exp(-((t-(r+0.035))/0.020).^2) ... % S wave
        + 0.30*exp(-((t-(r+0.25))/0.075).^2);     % T wave
end

baselineWander = 0.18*sin(2*pi*0.25*t);
powerlineNoise = 0.04*sin(2*pi*50*t);
whiteNoise = 0.06*randn(size(t));
noisyECG = cleanECG + baselineWander + powerlineNoise + whiteNoise;

% Baseline correction with a moving-average estimate; then smooth.
window = max(3, round(0.8*fs));
baselineEstimate = movmean(noisyECG, window);
baselineCorrected = noisyECG - baselineEstimate;
smoothWindow = 5;
filteredECG = movmean(baselineCorrected, smoothWindow);

% Toolbox-free local maximum detection with an amplitude threshold
% and a refractory period of 0.35 seconds.
threshold = 0.45 * max(filteredECG);
minDistance = round(0.35*fs);
candidate = find(filteredECG(2:end-1) > filteredECG(1:end-2) & ...
                 filteredECG(2:end-1) >= filteredECG(3:end) & ...
                 filteredECG(2:end-1) > threshold) + 1;
peakLocations = [];
for k = 1:numel(candidate)
    if isempty(peakLocations) || candidate(k)-peakLocations(end) >= minDistance
        peakLocations(end+1) = candidate(k); %#ok<SAGROW>
    elseif filteredECG(candidate(k)) > filteredECG(peakLocations(end))
        peakLocations(end) = candidate(k);
    end
end
peakTimes = t(peakLocations);
if numel(peakTimes) > 1
    RR = diff(peakTimes);
    heartRate = 60./RR;
    fprintf('Detected %d candidate R peaks. Mean rate: %.1f bpm\\n', ...
        numel(peakTimes), mean(heartRate));
else
    fprintf('Detected %d candidate R peaks. Adjust threshold if needed.\\n',numel(peakTimes));
end

figure('Color','w','Name','ECG noise removal and QRS detection');
subplot(3,1,1); plot(t,noisyECG); grid on; title('Noisy synthetic ECG');
ylabel('Amplitude (a.u.)');
subplot(3,1,2); plot(t,filteredECG); grid on; title('Baseline-corrected and smoothed ECG');
ylabel('Amplitude (a.u.)');
subplot(3,1,3); plot(t,filteredECG); hold on;
plot(peakTimes,filteredECG(peakLocations),'rv','MarkerFaceColor','r');
grid on; title('Detected candidate R peaks'); xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
outDir = fullfile(fileparts(mfilename('fullpath')),'..','results');
if ~exist(outDir,'dir'), mkdir(outDir); end
saveas(gcf,fullfile(outDir,'ecg_processing_result.png'));
