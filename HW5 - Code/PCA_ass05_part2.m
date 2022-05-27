clear all
close all
%% data
[Dp, dDp, ndist] = generate_sizehist; %initial distribution
gradp = [10^-12 1e-11 1e-10 1e-9]; %[atm]
N0 = 200*10^6; %[m^-3]
N = [50,200,500]*10^6; %[cm^-3]
gradp0 = 1e-9; %[atm]
%% Computation of Tau_{k}^{coag}
kappa = coag_loss_coef(N(3)*ndist,Dp);
tao_coag = 1./kappa/3600; %[h]
figure(1)
loglog(Dp,tao_coag,'.')
ylabel('\tau_{coag}')
xlabel('D_p')
xlim([10^-3 10^-1])
%% Computation of Tao_{k,k+1}^{cond}
gr = growth_rate(Dp*10^-6,gradp0);
tao_cond = -dDp(1:end-1)*10^-6./diff(gr)/3600;
figure(2)
loglog(Dp(1:end-1),tao_cond,'.')
ylabel('\tau_{cond}')
xlabel('D_p')
xlim([10^-3 10^-1])
ylim([10^-1 10^3])
%% Computation of Pr_{k->k+1}
Pr = exp(-tao_cond./tao_coag(1:end-1));
figure(3)
semilogx(Dp(1:end-1),Pr,'.')
xlim([10^-3 10^-1])


%% Variable Delta P
fig_gradP = figure(4);
for i = 1:length(gradp)
    subplot(length(gradp),2,2*i-1)
    hold on
    % TAO COAGULATION
        kappa = coag_loss_coef(N0*ndist,Dp);
        tao_coag = 1./kappa/3600; %[h]
        plot(Dp,tao_coag,'LineWidth',1.5)
    % TAO COAGULATION
        gr = growth_rate(Dp*10^-6,gradp(i));
        tao_cond = -dDp(1:end-1)*10^-6./diff(gr)/3600;
        plot(Dp(1:end-1),tao_cond,'LineWidth',1.5)
    hold off
    ylabel('Timescale [h]')
    xlabel('D_p [\mu m]')
    xlim([10^-3 10^-1]);
    ylim([10^-1 10^6]);
    yticks([10^-1 10^0 10^1 10^2 10^3 10^4 10^5])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend('\tau_{coag}','\tau_{cond}','location','Northwest')
    title('\nabla p_A = '+string(gradp(i)))
    % TAO CONDENSATION
    subplot(length(gradp),2,2*i)
    hold on
        Pr = exp(-tao_cond./tao_coag(1:end-1));
        plot(Dp(1:end-1),Pr,'LineWidth',1.5)
    xlim([10^-3 10^-1])  
    set(gca,'xscale','log');
    ylabel('Probability [-]')
    xlabel('D_p [\mu m]')
    title('\nabla p_A = '+string(gradp(i)))
    hold off  
end

%% Variable N
fig_N = figure(5);
for i = 1:length(N)
    subplot(length(N),2,2*i-1)
    hold on
    % TAO COAGULATION
        kappa = coag_loss_coef(N(i)*ndist,Dp);
        tao_coag = 1./kappa/3600; %[h]
        plot(Dp,tao_coag,'LineWidth',1.5)
    % TAO COAGULATION
        gr = growth_rate(Dp*10^-6,gradp0);
        tao_cond = -dDp(1:end-1)*10^-6./diff(gr)/3600;
        plot(Dp(1:end-1),tao_cond,'LineWidth',1.5)
    hold off
    ylabel('Timescale [h]')
    xlabel('D_p [\mu m]')
    xlim([10^-3 10^-1]);
    ylim([10^-1 10^6]);
    yticks([10^-1 10^0 10^1 10^2 10^3 10^4 10^5])
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    legend('\tau_{coag}','\tau_{cond}','location','Northwest')
    title('N = '+string(N(i)*10^-6))
    % TAO CONDENSATION
    subplot(length(N),2,2*i)
    hold on
        Pr = exp(-tao_cond./tao_coag(1:end-1));
        plot(Dp(1:end-1),Pr,'LineWidth',1.5)
    xlim([10^-3 10^-1])  
    set(gca,'xscale','log');
    ylabel('Probability [-]')
    xlabel('D_p [\mu m]')
    title('N = '+string(N(i)*10^-6))
    hold off  
end

%% SAVE FILES
% saveas(fig_gradP,'1.2_varGradP.png')
% saveas(fig_N,'1.2_varN.png')