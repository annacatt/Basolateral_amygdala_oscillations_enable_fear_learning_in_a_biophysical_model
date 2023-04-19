clear all

trialCOUNT = 1; 

b = [0, 0.4470, 0.7410];
o = [0.8500, 0.3250, 0.0980];
y = [0.9290, 0.6940, 0.1250];
p = [0.4940, 0.1840, 0.5560];
g = [0.4660, 0.6740, 0.1880];
a = [0.3010, 0.7450, 0.9330];
r = [0.6350, 0.0780, 0.1840];


%% DATA with rule biased towards depression (with no plastic F -> VIP)

external_stimulation = 3; %3: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
condition = 1; % 0: before FC; 1: during FC; 2: after FC;
%experiment = ...; %0: PV and SOM; 1: without SOM; 2: without PV; 3: neither PV nor SOM; 4: only PV+PN connections; 5: no connections at all; 6: no VIP; 7: with only PV; 8: with only SOM
stdp_rule = 1; %0: asym; 1: asym with larger tau and same amplitude; 2: sym; 3: no rule; 4:no rule but still evolution in time of weight as if there was plasticity

path = '';
filename = 'onlyPVlowexc';


%%
experiment_ = 7; % 7 for only PV (at low excitation level)
trial = 0;

for k = 1:length(experiment_)
    experiment = experiment_(k);

    M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);


    % For Fig. 4 main text
    %in = 6000/0.05;
    %en = 7000/0.05;

    % For Fig. S1C
    in = 6200/0.05;
    en = 6700/0.05;

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,4),'color',g,'LineWidth',1.5)
    hold on
    plot(M(in:en,1),M(in:en,6),'color',o,'LineWidth',1.5)
    legend('PV','F')
    xlim([in*0.05 en*0.05])

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,4),'color',g,'LineWidth',1.5)
    ylabel('PV')
    ylim([-100 50])
    set(gca, 'fontsize', 14)
    xlim([in*0.05 en*0.05])

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,5),'color',a,'LineWidth',1.5)
    hold on
    plot(M(in:en,1),M(in:en,6),'color',o,'LineWidth',1.5)
    ylim([-100 50])
    legend('ECS','F')
    set(gca, 'fontsize', 14)
    xlim([in*0.05 en*0.05])

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,9),'r')
    hold on
    plot(M(in:en,1),M(in:en,10),'b')
    legend('P: potentiation', 'M: depression')
    xlim([in*0.05 en*0.05])

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,9),'r')
    ylim([0 0.01])
    xlim([in*0.05 en*0.05])

    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),M(in:en,10),'b')
    ylim([-0.011 0])
    xlim([in*0.05 en*0.05])
    %legend('P: potentiation', 'M: depression')

    tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
    figure('Renderer', 'painters', 'Position', [10 10 500 100])
    plot(M(in:en,1),tmp(in:en))
    xlim([in*0.05 en*0.05])

    figure, plot(M(:,1),tmp(:,1))
end


%%

experiment_ = 1; % no SOM
trial = 0;

path = '';
filename = '';

M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment_),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);

% For Fig. 4 main text
in = 12200/0.05;%160000;
en = 13200/0.05;%170000;

% For Fig. S1D
in = 12700/0.05;%160000;
en = 13200/0.05;%170000;

figure,
plot(M(in:en,1),M(in:en,2))
set(gca, 'fontsize', 14)
xlim([in*0.05 en*0.05])

figure,
plot(M(in:en,1),M(in:en,5),'color',a,'LineWidth',1.5)
hold on
plot(M(in:en,1),M(in:en,6),'color',o,'LineWidth',1.5)
legend('ECS','F')
xlabel('[ms]')
set(gca, 'fontsize', 14)
xlim([in*0.05 en*0.05])

figure,
plot(M(in:en,1),M(in:en,10),'b')
ylim([-0.01 0])
xlabel('[ms]')
xlim([in*0.05 en*0.05])

figure,
plot(M(in:en,1),M(in:en,9),'r')
ylim([0 0.01])
xlabel('[ms]')
legend('potentiation')
set(gca, 'fontsize', 14)
xlim([in*0.05 en*0.05])


tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment_),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
figure, plot(M(in:en,1),tmp(in:en,1))
ylim([0 0.02])
xlim([in*0.05 en*0.05])


figure, plot(M(1:end,1),tmp(1:end,1))


