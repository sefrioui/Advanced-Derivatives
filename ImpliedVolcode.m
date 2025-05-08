
clear all
close all
clc

data = readtable("SX5E_Impliedvols.xlsx");
K = table2array(data(2:end,1))';
T = table2array(data(1,2:end));
iv = table2array(data(2:end,2:end));
t = 0;
q = 0;
r = 0;
T = [0, T];
S0 = 100;


K = K * S0; % The strike is scaled
dK = K(2) - K(1);  % Fixed strike spacing

% intialize call matrix 
C = zeros(length(K), length(T));
for i = 1:length(K)
    C(i, 1) = max(S0 - K(i), 0);  
end

% Call prices for each maturity and strike
for i = 1:length(K)
    for j = 1:length(T) - 1
        if iv(i, j) > 0
            C(i, j+1) = BSprice(S0, K(i), r, T(j+1), iv(i, j), q);
        end
    end
end

% Use A-H method to calibrate the vol
local_volatility = cell(1, length(T) - 1);
for j = 1:length(T) - 1
    opt_function = @(lv_params) sum((logical(iv(:, j)) .* ...
        (andreasen_price(T, K, j, C, lv_params, iv) - C(:, j+1))).^2);
    
    % Initial guess for local volatility parameters
    initial_guess = 0.3 * S0 * ones(1, length(K));
    % Optimize to get local volatility parameters
    [local_volatility{j}, ~] = fminsearch(opt_function, initial_guess);
    % Update call prices with calibrated local volatility
    C(:, j+1) = andreasen_price(T, K, j, C, local_volatility{j}, iv);
end

% Implied volatility matrix coputed 
implied_vol = zeros(size(iv));
for i = 1:length(K)
    for j = 1:length(T) - 1
        % Objective for implied volatility optimization
        vol_objective = @(vol) (BSprice(S0, K(i), r, T(j+1), vol, q) - C(i, j+1))^2;
        % Use fminsearch to estimate implied volatilities
        implied_vol(i, j) = fminsearch(vol_objective, 0.5);
    end
end

% Plot implied volatility surface
figure;
[x_mesh, y_mesh] = meshgrid(T(2:end), K);
surf(x_mesh, y_mesh, implied_vol);
xlabel('Time to Maturity (T)');
ylabel('Strike Price (K)');
zlabel('Implied Volatility');
title('Implied Volatility Surface');

% Interpolate for intermediate maturities T = 1 and T = 1.5
intermediate_maturities = [1, 1.5];
% Call interpolate_implied_vol with S0, r, and q
implied_vol_intermediate = interpolate_implied_vol(T, K, C, implied_vol, intermediate_maturities, S0, r, q);

% Plot final implied volatility surface with intermediate maturities
complete_vol_surface = [implied_vol, implied_vol_intermediate];
complete_maturities = [T(2:end), intermediate_maturities];
figure;
[x_full, y_full] = meshgrid(complete_maturities, K);
surf(x_full, y_full, complete_vol_surface);
xlabel('Time to Maturity (T)');
ylabel('Strike Price (K)');
zlabel('Implied Volatility');
title('Complete Implied Volatility Surface with Interpolated Maturities');

% Function Definitions

% Andreasen-Huge price computation function
function price_model = andreasen_price(T, K, j, current_prices, lv_params, iv)
    dT = T(j+1) - T(j);
    dK = K(2) - K(1);
   
    local_vol = zeros(length(K), 1);
    active_strikes = find(iv(1:end-1, j));
    idx = 1;
    for i = 1:length(K) - 1
        if i <= active_strikes(idx)
            local_vol(i) = lv_params(idx);
        else
            idx = min(idx + 1, length(active_strikes));
            local_vol(i) = lv_params(idx);
        end
    end

    % Compute z values for Dupire equation
    z = 0.5 * dT / dK^2 * local_vol.^2;
    
    % Assemble matrix A for solving Dupire equation
    A = diag([0, 2 * z(2:end-1)', 0]) + diag([-z(2:end-1)', 0], -1) + diag([0, -z(2:end-1)'], 1);
    A = A + eye(length(K));
    
    % Solve for new prices at expiration T(j+1)
    price_model = A \ current_prices(:, j);
end

function interpolated_vol = interpolate_implied_vol(T, K, current_prices, implied_vol, target_maturities, S0, r, q)
    % Function to interpolate implied volatilities for new intermediate maturities
    % using Andreasen-Huge method (second phase) as described in the slide.
    
    % Initialize the matrix for interpolated implied volatilities
    interpolated_vol = zeros(length(K), length(target_maturities));
    
    for k = 1:length(target_maturities)
        % Find the last maturity before the target maturity (i.e., T_i such that T_i <= target_maturities(k))
        idx = find(T <= target_maturities(k), 1, 'last');
        
        % Calculate the time step (Delta T) for interpolation
        dt = target_maturities(k) - T(idx);
        
        % Assume constant volatility between two market maturities
        % No need to recalibrate local volatility parameters
        z_values = 0.5 * dt / (K(2) - K(1))^2 * implied_vol(:, idx).^2;
        
        % Construct matrix A for solving prices at the intermediate maturity
        A = diag([0, 2 * z_values(2:end-1)', 0]) + diag([-z_values(2:end-1)', 0], -1) + diag([0, -z_values(2:end-1)'], 1);
        A = A + eye(length(K));
        
        % Solve for the option prices at the intermediate maturity using matrix A
        interpolated_prices = A \ current_prices(:, idx);
        
        % Since we assume constant volatility between two maturities, we do not recalibrate
        % Use the Black-Scholes formula directly to find implied volatility
        for i = 1:length(K)
            optim_fn = @(sigma) (BSprice(S0, K(i), r, target_maturities(k), sigma, q) - interpolated_prices(i))^2;
            interpolated_vol(i, k) = fminsearch(optim_fn, 0.5);
        end
    end
end

% Black-Scholes price calculation function
function price = BSprice(S, K, r, T, sigma, q)
    d1 = (log(S / K) + (r - q + 0.5 * sigma^2) * T) / (sigma * sqrt(T));
    d2 = d1 - sigma * sqrt(T);
    price = S * exp(-q * T) * normcdf(d1) - K * exp(-r * T) * normcdf(d2);
end


