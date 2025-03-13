% Clear all variables and close all plots
clear all
close all

% Create the 'images' folder if it does not exist
output_folder = 'fractional_images';
data_folder = 'fractional_data';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end

tic
% Define the time vector
t = linspace(0, 8, 500);  % Time vector

% Define the function f(t)

f = exp(-(t-3).^2)+0.025*(t+1);  % Example function


% Define the fractional derivative order values
alpha_values = [0.3 0.5 0.8 1];

% Initialize arrays to store the computed derivatives
caputo_derivatives = cell(length(alpha_values), 1);
rl_derivatives = cell(length(alpha_values), 1);
confor_derivatives = cell(length(alpha_values),1);
conori_derivatives = cell(length(alpha_values),1);

% Calculate the derivative of f(t) with respect to t
d_f = gradient(f, t);

% Calculate Caputo fractional derivatives
for i = 1:numel(alpha_values)
    alphav = alpha_values(i);
    delta = f(1) * (t-0).^(0-alphav) / gamma(0+1-alphav);
    delta(1) = ([0, 0, 1] * ([t(2:4).^2', t(2:4)', [1, 1, 1]'] \ (delta(2:4))'));
    
    caputo_derivatives{i} = numerical_caputo_derivative(d_f, t, alphav);
    rl_derivatives{i} = caputo_derivatives{i} + delta;
    confor_derivatives{i} = exp(t*(1-alphav)) .* d_f;
    conori_derivatives{i} = t.^(1-alphav) .* d_f;

    % Save data to .dat files
    save_data([t', caputo_derivatives{i}'], fullfile(data_folder, ['caputoalpha' num2str(alphav) '.dat']));
    save_data([t', rl_derivatives{i}'], fullfile(data_folder, ['rlalpha' num2str(alphav) '.dat']));
    save_data([t', confor_derivatives{i}'], fullfile(data_folder, ['conforalpha' num2str(alphav) '.dat']));
    save_data([t', conori_derivatives{i}'], fullfile(data_folder, ['conorialpha' num2str(alphav) '.dat']));
end

% Plot all fractional derivatives
figure;
hold on;
plot(t, f, 'k', 'LineWidth', 2, 'DisplayName', 'Original Function');
plot(t, caputo_derivatives{2}, 'LineWidth', 1.5, 'DisplayName', ['Caputo, \alpha=' num2str(alpha_values(2))]);
plot(t, rl_derivatives{2}, '--', 'LineWidth', 1.5, 'DisplayName', ['Riemann-Liouville, \alpha=' num2str(alpha_values(2))]);
plot(t, confor_derivatives{2}, '-.', 'LineWidth', 1.5, 'DisplayName', ['Conformable, \alpha=' num2str(alpha_values(2))]);
plot(t, conori_derivatives{2}, '-.', 'LineWidth', 1.5, 'DisplayName', ['Orig. Conformable, \alpha=' num2str(alpha_values(2))]);
xlabel('t');
ylabel('Fractional Derivatives');
legend;
title('Comparison of Fractional Derivatives');
hold off;

function save_data(data, filename)
    fid = fopen(filename, 'w');
    fprintf(fid, '%.6f %.6f\n', data');
    fclose(fid);
end

% Define a function to compute the numerical Caputo fractional derivative
function CD = numerical_caputo_derivative(f, t, alpha)
    n = ceil(alpha); % Round alpha up to the nearest integer
    CD = zeros(size(f)); % Initialize the array to store the result

    % Loop over each time point
    for j = 1:numel(t)
        integ = 0; % Initialize the integral value

        % Compute the integral using summation
        for k = 1:j
            if k == 1
                f_diff = f(1); % Compute the difference at the first point
            else
                f_diff = f(k) - f(k-1); % Compute the difference at other points
            end
            % Compute the integrand and add it to the integral value
            integ = integ + ((t(j) - t(k))^(n-alpha)) * f_diff / gamma(n-alpha+1);
        end

        % Store the computed Caputo fractional derivative value
        CD(j) = integ;
    end
end
