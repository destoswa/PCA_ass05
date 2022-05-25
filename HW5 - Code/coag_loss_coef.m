%% -----------------------------------------------------------------------------

 function kappa = coag_loss_coef(n, Dp)

 kmax = length(Dp);
 kappa = zeros(size(Dp));

 for k = 1:kmax
     j = (k+1):kmax;
     for i = 1:j
        kappa(k) = kappa(k) + n(k)*coag_kernel(Dp(k),Dp(i))*n(i);
     end
 end

 %% -----------------------------------------------------------------------------