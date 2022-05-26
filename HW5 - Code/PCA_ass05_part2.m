clear all
close all
%% data
[Dp, dDp, ndist] = generate_sizehist;
gradp = [1e-12 1e-11 1e-10 1e-9]; %[atm]
N0 = 200; %[cm^-3]
N = [50,200,500]; %[cm^-3]
gradp0 = 1e-9; %[atm]
%% Computation of Tau_{k}^{coag}
tao_coag = 1./coag_loss_coef(ndist,Dp);
figure(1)
semilogy(Dp,tao_coag,'.')
%% Computation of Tao_{k,k+1}^{cond}
tao_cond = zeros(length(Dp)-1,1);
for i=1:length(tao_cond)
    tao_cond(i) = Dp(i+1)-Dp(i)/(growth_rate(Dp(i+1),dDp(i+1))-growth_rate(Dp(i),dDp(i)));
end
figure(2)
semilogy(Dp(1:end-1),tao_cond,'.')
%% Computation of Pr_{k->k+1}