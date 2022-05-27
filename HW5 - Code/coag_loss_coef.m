%% -----------------------------------------------------------------------------

 function kappa = coag_loss_coef(n, Dp)

 kmax = length(Dp);
 kappa = zeros(size(Dp));
 Dp = Dp*10^-6; %[m]
 for k = 1:kmax
     j = (k+1):kmax;
     kappa(k) = 0.5*coag_kernel(Dp(k),Dp(k))*n(k);
     for i = (k+1):kmax
        kappa(k) = kappa(k) + coag_kernel(Dp(k),Dp(i))*n(i);
     end
 end

 %% -----------------------------------------------------------------------------