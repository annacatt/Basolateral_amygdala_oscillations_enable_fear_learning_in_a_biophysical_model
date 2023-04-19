clear all

trialCOUNT = 40; 

%% DATA with rule biased towards depression (with no plastic F -> VIP)

external_stimulation = 3; %0: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
condition = 1; % 0: before FC; 1: during FC; 2: after FC;
%experiment = 0; %0: PV and SOM; 1: without SOM; 2: without PV; 3: neither PV nor SOM; 4: only PV+PN connections; 5: no connections at all; 6: no VIP; 7: with only PV; 8: with only SOM
stdp_rule = 1; %0: asym; 1: asym with larger tau and same amplitude; 2: sym; 3: no rule; 4:no rule but still evolution in time of weight as if there was plasticity

path = '';
filename = 'singleneuronnet';

PLOT = 0;

experiment_ = 7%[0 6 1 2 7];% 3];%0%[0 6 1 2 3];%1; % without SOM;%0;%[0 6 1 2 3]; % 7 for only PV

for k = 1:length(experiment_)
    experiment = experiment_(k);
    for trial = 0:trialCOUNT-1
        [k trial]


        %%M = readmatrix(['/Users/ac/Documents/data/hodgkinhuxley_c++/STDP2/',filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);
        %M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);

        if PLOT == 1

            
            in = 7700/0.05; %7300/0.05;
            en = 8700/0.05; %8300/0.05;

            in = 8150/0.05; %7300/0.05;
            en = 8650/0.05; %8300/0.05;

            figure('Renderer', 'painters', 'Position', [10 10 500 100])
            plot(M(in:en,1),M(in:en,4),'color',g,'LineWidth',1.5)
            hold on
            plot(M(in:en,1),M(in:en,6),'color',o,'LineWidth',1.5)
            legend('PV','F')
            xlim([in*0.05 en*0.05])


            figure('Renderer', 'painters', 'Position', [10 10 500 100])
            plot(M(in:en,1),M(in:en,2),'color',b,'LineWidth',1.5)
            %hold on
            %plot(M(1:40000,1),M(1:40000,6)./10000)
            ylim([-100 50])
            ylabel('VIP')
            set(gca, 'fontsize', 14)
            xlim([in*0.05 en*0.05])

            figure('Renderer', 'painters', 'Position', [10 10 500 100])
            plot(M(in:en,1),M(in:en,3),'color',p,'LineWidth',1.5)
            ylabel('SOM')
            ylim([-100 50])
            set(gca, 'fontsize', 14)
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
            xlim([in*0.05 en*0.05])
            ylim([0 0.01])


            figure('Renderer', 'painters', 'Position', [10 10 500 100])
            plot(M(in:en,1),M(in:en,10),'b')
            %legend('P: potentiation', 'M: depression')
            xlim([in*0.05 en*0.05])
            ylim([-0.01 0])

            tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
            figure('Renderer', 'painters', 'Position', [10 10 500 100])
            plot(M(in:en,1),tmp(in:en))
            xlim([in*0.05 en*0.05])

        end


        tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);

        if experiment == 3 && trial == 2
            gAMPA(k).v(trial+1,:,:) = zeros(200001,1);
        else
            gAMPA(k).v(trial+1,:,:) = tmp;
        end

        if PLOT == 1
            figure, plot(M(:,1),tmp)
        end

    end

end

M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);

figure,
[lineOut(1), fillOut(1)] = stdshade(squeeze(gAMPA(1).v(:,:,1)),0.2,'r',M(:,1));%,alpha,acolor,F,smth);
hold on
[lineOut(2), fillOut(2)] = stdshade(squeeze(gAMPA(2).v(:,:,1)),0.2,'g',M(:,1));%,alpha,acolor,F,smth);
[lineOut(3), fillOut(3)] = stdshade(squeeze(gAMPA(3).v(:,:,1)),0.2,'b',M(:,1));%,alpha,acolor,F,smth);
[lineOut(4), fillOut(4)] = stdshade(squeeze(gAMPA(4).v(:,:,1)),0.2,'k',M(:,1));%,alpha,acolor,F,smth);ù[lineOut(4), fillOut(4)] = stdshade(gAMPA(4).v,0.2,'g',M(:,1));%,alpha,acolor,F,smth);ù
%[lineOut(5), fillOut(5)] = stdshade(gAMPA(5).v,0.2,'m',M(:,1));%,alpha,acolor,F,smth);
hold on
legend(lineOut,'whole network','no VIP', 'no SOM', 'no PV')%, 'no SOM and PV')
xlabel('t [ms]')
ylabel('w change')
set(gca, 'fontsize', 14)


