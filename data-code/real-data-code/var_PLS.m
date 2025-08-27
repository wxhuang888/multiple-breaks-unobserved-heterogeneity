function [vari] = var_PLS(T, a, y, X)

% Number of cross-sectional units
Nk = size(y, 1) / T;

% Calculate Phi (X'X scaled by Nk * T)
Phi = 1 / (Nk * T) * (X' * X);

% Residuals
ui = y - X * a;

% Variance of residuals (assuming homoskedasticity)
sigma2 = (ui' * ui) / (Nk * T - size(X, 2)); % Divided by degrees of freedom

% Simplified variance-covariance matrix
vari = sigma2 * (1 / (Nk * T)) * pinv(Phi);

end