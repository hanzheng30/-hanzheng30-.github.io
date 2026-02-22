% Elliptic Robot Parameters
m = 1; g = 9.81;
ae = 0.15; be = 0.05; d = 0.4;

% Height and horizontal offset of CoM relative to contact point
h = @(theta) sqrt(ae^2 * sin(theta).^2 + be^2 * cos(theta).^2);
xcom = @(theta) ((ae^2 - be^2) * sin(theta) .* cos(theta)) ./ h(theta);

% Arc length (rolling distance on floor)
xs = @(theta) arrayfun(@(t) integral(@(phi) ...
    sqrt(ae^2 * cos(phi).^2 + be^2 * sin(phi).^2), 0, t), theta);

% Tail length L from tip (-d, 0) to CoM (xs - xcom, h)
L = @(theta) sqrt((xs(theta) - xcom(theta) + d).^2 + h(theta).^2);

% Lever arm of Force F about contact point
% Lever = |(CoM - Contact) x (Tail Vector)| / L
lever = @(theta) (h(theta) .* (xs(theta) + d)) ./ L(theta);

% Force and Work Integrand
F = @(theta) (m * g * xcom(theta)) ./ lever(theta);
delta = 1e-6;
dL_dtheta = @(theta) (L(theta + delta) - L(theta)) / delta;
work_integrand = @(theta) F(theta) .* dL_dtheta(theta);

% Calculate Work to tipping point (pi/2)
W_calc = integral(work_integrand, 0, pi/2);
fprintf('Updated Calculated Work: %.3f J\n', W_calc);