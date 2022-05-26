%import coag_kernel.m
%import generate_sizehist.m
%import growth_rate.m
close all
clear all
%% integrate over {0, 1, 2, ..., 18000} seconds with an initial diameter of 3 nm 
%%   for pressure gradients of 1e-12 and 1e-10 atmospheres.

%% define values
t = linspace(0, 5) .* 3600; % [s]
Dp0 = 3e-9; % [m]
gradP1 = 1e-12; % [atm] equivalent to 1 ppt
gradP2 = 1e-10; % [atm] equivalent to 100 ppt

if(exist('OCTAVE_VERSION', 'builtin') ~= 0)
  %% integrate (Octave)
  y1 = lsode(@(Dp, t) growth_rate(Dp, gradP1), Dp0, t);
  y2 = lsode(@(Dp, t) growth_rate(Dp, gradP2), Dp0, t);
else
  %% integrate (MATLAB)
  [~,y1] = ode45(@(t, Dp) growth_rate(Dp, gradP1), t, Dp0);
  [~,y2] = ode45(@(t, Dp) growth_rate(Dp, gradP2), t, Dp0);
end

%% plot
plot(t / 3600, [y1,y2] * 1e6)
set(gca, 'yscale', 'log')
xlabel('Time [h]')
ylabel('D_p [\mum]')
legend({'\nabla p_A = 10^{-12} atm', '\nabla p_A = 10^{-10} atm'}, 'location', 'northwest')