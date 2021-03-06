clear all
close all
%% data
[Dp, dDp, ndist] = generate_sizehist; %initial distribution
gradp = [1e-12 1e-11 1e-10 1e-9]; %[atm]
N0 = 200*10^6; %[m^-3]
N = [50,200,500]; %[cm^-3]
gradp0 = 1e-12; %[atm]
%% Computation of Tau_{k}^{coag}
kappa = coag_loss_coef(N0*ndist,Dp);
tao_coag = 1./kappa/3600; %[h]
figure(1)
loglog(Dp,tao_coag,'.')
ylabel('\tau_{coag}')
xlabel('D_p')
xlim([10^-3 10^-1])
%% Computation of Tao_{k,k+1}^{cond}
% tao_cond = zeros(length(Dp)-1,1);
gr = growth_rate(Dp*10^-6,gradp0);
tao_cond = -dDp(1:end-1)*10^-6./diff(gr)/3600;
% for i=1:length(tao_cond)
%     tao_cond(i) = -((Dp(i+1)-Dp(i))*10^-6)/(growth_rate(Dp(i+1)*10^-6,gradp0)-growth_rate(Dp(i)*10^-6,gradp0))/3600; % [h]
% end
figure(2)
% loglog(Dp(1:end-1),tao_cond,'.')
loglog(Dp(1:end-1),tao_cond,'.')
ylabel('\tau_{cond}')
xlabel('D_p')
xlim([10^-3 10^-1])
ylim([10^0 10^3])
%% Computation of Pr_{k->k+1}


%% Computation of Pr_{3->100nm}
Pr=zeros(199,1);
Pr(1)=exp(-tao_cond(1)/tao_coag(1));
for k=1:198 
   Pr(k+1)=Pr(k)*exp(-tao_cond(k+1)/tao_coag(k+1));
end

figure(3)
% loglog(Dp(1:end-1),tao_cond,'.')
loglog(Dp(1:end-1),Pr,'.')
ylabel('P_{r}')
xlabel('D_p')
xlim([10^-3 10^-1])
ylim([1 0])
