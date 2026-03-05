% 1. Setup time vector
fs = 500;                    % Sampling frequency (Hz)
t = 0:1/fs:1;                % 1 second duration

% 2. Create the modulated signal (10Hz sine wave fading in)
carrier = sin(2 * pi * 10 * t);
modulation = t;              % The "envelope" we want to extract
signal = modulation .* carrier;

% 3. Extract the envelope using the Hilbert Transform
analytic_signal = hilbert(signal);    % Returns s(t) + j*H{s(t)}
envelope = abs(analytic_signal);      % Magnitude of the complex signal

% 4. Visualize the results
figure;
plot(t, signal, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Original Signal');
hold on;
plot(t, envelope, 'r', 'LineWidth', 2, 'DisplayName', 'Hilbert Envelope');
plot(t, abs(signal), 'g--', 'DisplayName', 'Absolute Value');

title('Hilbert Envelope vs. Simple Absolute Value');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Location', 'best');
grid on;
