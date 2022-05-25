%% -----------------------------------------------------------------------------

 function rate = growth_rate(D_p, gradP)

 %% define constants
M = 98.079; %[g/mol];
D_AB = 0.1*10^-4; %[m^2/s];
rho = 1.0*10^6; %[g/(m^3)];
R = 8.20573660809596*10^-5; %[m^3*atm/(K*mol)];
T = 298; %[K];
p_A = 1; %[atm];
alpha = 1;
lambda = 0.118*10^-6; %[m];
Kn = 2*lambda/D_p; 

 %% calculate Fuchs-Sutugin correction factor
 f = (0.75*alpha*(1+Kn))/(Kn^2+Kn+0.283*Kn*alpha+0.75*alpha);

 %% calculate I(Dp)
 rate = (4*D_AB*M*f*gradP)/(R*T*rho*D_p); %[m/s]

 end

 %% -----------------------------------------------------------------------------