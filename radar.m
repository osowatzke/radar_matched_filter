classdef radar < handle

    properties (Constant)
        % Speed of Light (m/s)
        c = 3e8;

        % Boltzmann's Constant (J/K)
        k = 1.38e-23;
    end

    properties
        % Carrier Frequency (Hz)
        fc = 10e9;

        % Sample Rate (Hz)
        fs = 1e8;

        % Transmit Power (dBW)
        Pt = 20;

        % Antenna Gain (dB)
        G = 10*log10(26000);

        % Temperature (K)
        T = 290;

        % Noise Figure (dB)
        F = 10;

        % System Loss (dB)
        Ls = 5;

        % duty cycle
        duty_cyle = 0.2;

        % Pulse Repetition Frequency (Hz)
        PRF = 50e3;

        % Bandwidth of chirp (Hz)
        B = 1e8;
    end

    methods

        % class constructor allows for parameters to be set using comma
        % separated list. Ex: r = radar('fc', 1e9) creates a radar class 
        % with fc (Carrier Frequency) property set to 1e9.
        function self = radar(varargin)
            for i = 1:2:nargin
                self.(varargin{i}) = varargin{i+1};
            end
        end

        % function computes PRI in s
        function y = PRI(self)
            y = 1/self.PRF;
        end

        % function computes wavelength in m
        function y = lambda(self)
            y = self.c/self.fc;
        end

        % function computes length of chirp pulse in s
        function y = tau(self)
            y = self.PRI*self.duty_cyle;
        end

        % function computes sample period of ADC in s
        function y = Ts(self)
            y = 1/self.fs;
        end

        % Compute time delay of target in s given range in m
        function td = compute_time_delay(self, R)
            td = 2*R/self.c;
        end

        % Compute range of target in m given range in s
        function R = compute_range(self, td)
            R = self.c*td/2;
        end

        % compute target radar cross section required to generate signal 
        % with desired SNR
        function sigma = compute_RCS(self, R, SNR)

            % compute SNR corresponding to RCS of 1 m^2
            SNR0 = self.compute_target_SNR(R, 1);

            % compute required RCS
            sigma = 10.^((SNR - SNR0)/10);
        end    

        % function performs generates ADC signal
        % and runs signal through matched filter
        function [adc_out, adc_SNR, mf_out, mf_SNR, mf_SNRt] = ...
            run_scenario(self, R, sigma)

            % compute ADC SNR
            adc_SNR = self.compute_target_SNR(R, sigma);

            % compute ADC signal
            [adc_out, x, n] = self.generate_received_signal(R, adc_SNR);

            % take maximum SNR if there are multiple targets
            adc_SNR = max(adc_SNR);

            % run matched filter
            [mf_out, x, n] = self.generate_matched_filter_output(x, n);

            % compute matched filter SNR
            mf_SNR = self.compute_mf_SNR(x, n);

            % compute theoretical matched filter SNR
            mf_SNRt = self.compute_theoretical_mf_SNR(adc_SNR);
        end
    end

    methods (Access=protected)

        % function computes transmit code
        function y = code(self)

            % Generate array of transmit times
            t = 0:self.Ts:(self.tau-self.Ts);

            % Generate Transmit Pulse
            y = exp(1i*pi*self.B*t.^2/self.tau);
        end

        % function computes matched filter
        function y = matched_filter(self)

            % matched filter is flip and conjugate of code
            y = flip(conj(self.code));
        end

        % compute atmospheric loss (dB)
        function y = La(self, R)

            % atmospheric pressure (Pa)
            P = 101325;

            % water vapor density (g/m^3)
            wvd = 7.5;

            % compute atmospheric loss (dB)
            y = gaspl(R, self.fc, self.T, P, wvd);
        end

        % function computes target SNR
        function SNR = compute_target_SNR(self, R, sigma)
            
            SNR = ...
                self.Pt + ... % Transmit Power (dB)
                2*self.G + ... % Antenna Gain (dB)
                20*log10(self.lambda) + ... % lambda^2 (dB)
                10*log10(sigma(:)) - ... % Radar Cross Section (dB)
                30*log10(4*pi) - ... % (4*pi)^3 (dB)
                40*log10(R(:)) - ... % R^4 (dB)
                10*log10(self.k) - ... % Boltzmann's Constant (dB)
                10*log10(self.T) - ... % Temperature (dB)
                10*log10(self.fs) - ... % Bandwidth (dB)
                self.F - ... % Noise Figure (dB)
                self.Ls - ... % System Loss (dB)
                self.La(R(:)); % Atmospheric Loss (dB)
        end

        % function computes transmit signal
        function y = generate_transmit_signal(self)

            % function computes length of transmit signal in samples
            len_tx_sig = round(self.PRI * self.fs);

            % transmit signal is code plus listening period
            y = [self.code zeros(1,len_tx_sig - length(self.code))];
        end

        % function generates received ADC signal
        function [y, x, n] = generate_received_signal(self, R, SNR)

            % generate array of times
            t = 0:self.Ts:(self.PRI-self.Ts);

            % generate complex noise of unit power
            n = 1/sqrt(2)*randn(1,length(t)) + 1i*1/sqrt(2)*randn(1,length(t));

            % compute delay of each target return (s)
            td = self.compute_time_delay(R);

            % generate transmit signal
            tx_sig = self.generate_transmit_signal;

            % convert time delay to samples
            % cannot delay more than length of signal
            td = round(td/self.Ts);

            % empty array for received signal
            x = zeros(1,length(t));

            % add each target return
            for i = 1:length(td)

                % determine how much signal needs to be delayed
                % cannot delay more than length of signal
                delay_samp = min([td(i) length(tx_sig)]);

                % compute unit power received signal for a given target
                rx_sig = [zeros(1,delay_samp) tx_sig(1:(end-delay_samp))];

                % scale target return to produce desired SNR
                rx_sig = rx_sig * (10^(SNR(i)/20));

                % add to return signal
                x = x + rx_sig;
            end

            % compute length of blanking
            len_code = length(self.code);

            % Let user know if any of there target returns are eclipsed
            if sum(td < len_code) > 0
                warning('One or More Targets is Partially Eclipsed');
            end

            % Add blanking
            n(1:len_code) = 0;
            x(1:len_code) = 0;

            % compute return signal (target returns + noise)
            y = x + n;
        end

        % function generates matched filter output
        function [y, x, n] = generate_matched_filter_output(self, x, n)

            % run signal through matched filter
            x = filter(self.matched_filter,1,x);

            % run noise through matched filter
            n = filter(self.matched_filter,1,n);

            % matched filter output is signal + noise
            y = x + n;
        end

        % function computes SNR of matched filter output
        function SNR = compute_mf_SNR(self, x, n)

            % get Signal Power at peak of matched filter output
            Ps = max(abs(x).^2);

            % Number of Samples of chargeup
            % Blanking + Chargeup of Filter
            Nchargeup = 2*length(self.code) - 1;

            % get Noise Power after matched filter
            Pn = mean(abs(n((Nchargeup+1):end).^2));

            % Compute SNR of matched filter output
            SNR = 10*log10(Ps/Pn);
        end

        % compute theoretical matched filter SNR
        function SNR = compute_theoretical_mf_SNR(self, adc_SNR)

            % Matched Filter Increases SNR by filter length
            SNR = adc_SNR + 10*log10(length(self.matched_filter));
        end
    end
end