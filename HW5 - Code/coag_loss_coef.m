%% -----------------------------------------------------------------------------

 function kappa = coag_loss_coef(n, Dp)

 kmax = length(Dp);
 kappa = zeros(size(Dp));

 for k = 1:kmax

     j = (k+1):kmax
     kappa(k) = {... you will need to call coag_kernel() ...}

 end

 %% -----------------------------------------------------------------------------