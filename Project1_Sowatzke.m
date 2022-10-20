% Author: O. Sowatzke
%
% Date: 10/20/2022
%
% Subject: Script finds the SNR of a single target return before and after
% matched filtering. Script also examines the case when two targets are
% resolvable and when two targets are not resolvable
%

%% Preprocessing

% generate a radar object
r = radar;

% generate array of times for plotting 
t = 0:r.Ts:(r.PRI-r.Ts);

%% Find SNR before and after Matched Filter

% compute time delay (s) corresponding to 25% of PRI
td = 0.25*r.PRI;

% compute range (m) required to place target at distance corresponding to 
% time delay
R = r.compute_range(td);

% Radar cross section of target (m^2) to provide a 20 dB of SNR at ADC
sigma = r.compute_RCS(R,20);

% Output Configuration of Scenario
fprintf('\nScenario 1:\n\n')
fprintf('\tTarget 1:\n');
fprintf('\t\tR = %g m\n', R);
fprintf('\t\tsigma = %g m^2\n', sigma);

% Generate Targets and Run Matched Filter
[adc_sig, adc_SNR, mf_sig, mf_SNR, mf_SNRt] = r.run_scenario(R, sigma);

% plot ADC Signal
adc_sig_norm = adc_sig./max(abs(adc_sig));
figure(1);
clf;
plot(t, db(adc_sig_norm));
xlabel('Time (s)');
ylabel('Magnitude (dB)');
title(sprintf('ADC Signal : SNR = %.2f dB', adc_SNR));
xlim([t(1) t(end)]);
ylim([-80 0]);
grid on;

% Plot Matched Filter Output
mf_sig_norm = mf_sig./max(abs(mf_sig));
figure(2);
clf;
plot(t, db(mf_sig_norm));
xlabel('Time (s)');
ylabel('Magnitude (dB)');
title(sprintf('Matched Filter Output : SNR = %.2f dB : SNRt = %.2f dB', mf_SNR, mf_SNRt));
xlim([t(1) t(end)]);
ylim([-80 0]);
grid on;

%% Resolvable Targets (Large Separation)

% compute time delays (s) corresponding to 25% of PRI
% and 25% of PRI + 2 Code Lengths
td = 0.25*r.PRI;
td = [td, td + 2*r.duty_cyle*r.PRI];

% compute range (m) required to place target at distance corresponding to 
% time delays
R = r.compute_range(td);

% Radar cross section of targets (m^2) to provide a 20 dB of SNR at ADC
sigma = r.compute_RCS(R,20);

% Output Configuration of Scenario
fprintf('\nScenario 2:\n\n')
fprintf('\tTarget 1:\n');
fprintf('\t\tR = %g m\n', R(1));
fprintf('\t\tsigma = %g m^2\n', sigma(1));
fprintf('\tTarget 2:\n');
fprintf('\t\tR = %g m\n', R(2));
fprintf('\t\tsigma = %g m^2\n', sigma(2));

% Generate Targets and Run Matched Filter
[~, ~, mf_sig] = r.run_scenario(R, sigma);

% generate array of times for plotting 
t = 0:r.Ts:(r.PRI-r.Ts);

% Plot Matched Filter Output
mf_sig_norm = mf_sig./max(abs(mf_sig));
figure(3);
clf;
plot(t, db(mf_sig_norm));
xlabel('Time (s)');
ylabel('Magnitude (dB)');
title('Matched Filter Output (Resolvable Targets)');
xlim([t(1) t(end)]);
ylim([-80 0]);
grid on;

%% Resolvable Targets (Small Separation)

% compute time delays (s) corresponding to 25% of PRI
% and 25% of PRI + 2 ADC Samples
td = 0.25*r.PRI;
td = [td td+2*r.Ts];

% compute range (m) required to place target at distance corresponding to 
% time delays
R = r.compute_range(td);

% Radar cross section of targets (m^2) to provide a 20 dB of SNR at ADC
sigma = r.compute_RCS(R,20);

% Output Configuration of Scenario
fprintf('\nScenario 3:\n\n')
fprintf('\tTarget 1:\n');
fprintf('\t\tR = %g m\n', R(1));
fprintf('\t\tsigma = %g m^2\n', sigma(1));
fprintf('\tTarget 2:\n');
fprintf('\t\tR = %g m\n', R(2));
fprintf('\t\tsigma = %g m^2\n', sigma(2));

% Generate Targets and Run Matched Filter
[~, ~, mf_sig] = r.run_scenario(R, sigma);

% generate array of times for plotting 
t = 0:r.Ts:(r.PRI-r.Ts);

% Plot Matched Filter Output
mf_sig_norm = mf_sig./max(abs(mf_sig));
figure(4);
clf;
plot(t, db(mf_sig_norm));
xlabel('Time (s)');
ylabel('Magnitude (dB)');
title('Matched Filter Output (Resolvable Targets)');
xlim([t(1) t(end)]);
ylim([-80 0]);
grid on;

%% Unresolvable Targets
% compute time delays (s) corresponding to 25% of PRI
% and 25% of PRI + 1 ADC Sample
td = 0.25*r.PRI;
td = [td td+r.Ts];

% compute range (m) required to place target at distance corresponding to 
% time delays
R = r.compute_range(td);

% Radar cross section of targets (m^2) to provide a 20 dB of SNR at ADC
sigma = r.compute_RCS(R,20);

% Output Configuration of Scenario
fprintf('\nScenario 4:\n\n')
fprintf('\tTarget 1:\n');
fprintf('\t\tR = %g m\n', R(1));
fprintf('\t\tsigma = %g m^2\n', sigma(1));
fprintf('\tTarget 2:\n');
fprintf('\t\tR = %g m\n', R(2));
fprintf('\t\tsigma = %g m^2\n', sigma(2));

% Generate Targets and Run Matched Filter
[~, ~, mf_sig] = r.run_scenario(R, sigma);

% generate array of times for plotting 
t = 0:r.Ts:(r.PRI-r.Ts);

% Plot Matched Filter Output
mf_sig_norm = mf_sig./max(abs(mf_sig));
figure(5);
clf;
plot(t, db(mf_sig_norm));
xlabel('Time (s)');
ylabel('Magnitude (dB)');
title('Matched Filter Output (Unresolvable Targets)');
xlim([t(1) t(end)]);
ylim([-80 0]);
grid on;