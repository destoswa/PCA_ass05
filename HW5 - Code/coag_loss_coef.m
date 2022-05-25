%% -----------------------------------------------------------------------------

 function kappa = coag_loss_coef(n, Dp)

 kmax = length(Dp);
 kappa = zeros(size(Dp));

 for k = 1:kmax

     j = (k+1):kmax;
     kappa(k) = coag_kernel(Dp(k),Dp(j))*n(j);

 end

 %% -----------------------------------------------------------------------------