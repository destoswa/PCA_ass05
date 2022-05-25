function K_12 = coag_kernel(Dp_1, Dp_2)
  %% Dp_1, Dp_2 in [m]

  %% --- constants / coefficients ---
  kB = 1.3806e-23;	% m^2 kg s^-2 K^-1
  T = 298;		% K
  mu_B = 1.83e-5;	% kg m^-1 s^-1
  rho_p = 1e-3;		% kg m^-3
  lambda = 0.0686e-6;   % m, Seinfeld and Pandis 2006, Ch. 9

  %% --- diffusion coefficient of particle ---
  Kn = @(Dp) 2*lambda./Dp;
  [A, B, C] = deal(0.864, 0.29, 1.25);		%% Jacobson 2005, Chapter 15
  Cc = @(Dp) 1 + Kn(Dp).*(A + B.*exp(-C./Kn(Dp)));
  D = @(Dp) kB.*T.*Cc(Dp)./(3.*pi.*mu_B.*Dp);	% diffusion coefficient

  %% --- mean speed of particle ---
  mass = @(Dp) rho_p.*(pi/6).*Dp.^3;			% mass of spherical particle
  cbar_p = @(Dp) (8.*kB.*T./(pi.*mass(Dp))).^(1/2);	% mean thermal speed of particle in air

  %% --- correction factor ---
  elle = @(Dp) 8.*D(Dp)./(pi.*cbar_p(Dp));
  g = @(Dp) 1./(3.*Dp.*elle(Dp)) .* ((Dp+elle(Dp)).^3 - (Dp.^2+elle(Dp).^2).^(3/2)) - Dp;
  mix_12 = @(fn, Dp_1, Dp_2) (fn(Dp_1).^2 + fn(Dp_2).^2).^(1/2);
  beta = @(Dp_1, Dp_2) ((Dp_1 + Dp_2)./(Dp_1 + Dp_2 + 2.*mix_12(g, Dp_1, Dp_2)) + ...
			8.*(D(Dp_1) + D(Dp_2)) ./ (mix_12(cbar_p, Dp_1, Dp_2).*(Dp_1 + Dp_2))).^(-1);

  %% --- coagulation coefficient ---
  K_12 = 2*pi .* (D(Dp_1) + D(Dp_2)) .* (Dp_1 + Dp_2) .* beta(Dp_1, Dp_2); % m^3 s^-1

end