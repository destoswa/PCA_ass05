clear all
close all
%% data
[Dp, dDp, ndist] = generate_sizehist;
gradp = [1e-12 1e-11 1e-10 1e-9]; %[atm]
N0 = 200; %[cm^-3]
N = [50,200,500]; %[cm^-3]
gradp0 = 1e-9; %[atm]
%% Computation of Tau_{k}^{coag}
tao_coag = ndist./coag_loss_coef(ndist,Dp);
figure(1)
semilogy(Dp,tao_coag,'.')
%% Computation of Tao_{k,k+1}^{cond}

%% Computation of Pr_{k->k+1}