clear all

trialCOUNT = 1;

b = [0, 0.4470, 0.7410];
o = [0.8500, 0.3250, 0.0980];
y = [0.9290, 0.6940, 0.1250];
p = [0.4940, 0.1840, 0.5560];
g = [0.4660, 0.6740, 0.1880];
a = [0.3010, 0.7450, 0.9330];
r = [0.6350, 0.0780, 0.1840];

external_stimulation = 3; %0: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
condition = 1; % 0: before FC; 1: during FC; 2: after FC;
experiment = 0; %0: PV and SOM; 1: without SOM; 2: without PV; 3: neither PV nor SOM; 4: only PV+PN connections; 5: no connections at all; 6: no VIP; 7: with only PV; 8: with only SOM
stdp_rule = 1; %0: asym; 1: asym with larger tau and same amplitude; 2: sym; 3: no rule; 4:no rule but still evolution in time of weight as if there was plasticity

path = ''; %'/Volumes/TOSHIBA/fearconditioning/heterogeneous/';
filename = 'het';


PLOT = 1;

%% Fig. 5BC

experiment_ = 0; % full network
for k = 1:length(experiment_)
    experiment = experiment_(k);
    for trial = 0:trialCOUNT-1 %4
        [k trial]

        if PLOT == 1
            trial = 16;

            M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);

            t_length = 100000; % 5 seconds

            %VIP
            act_VIP = M(1:t_length,2:4);
            spike_VIP = cell(3,1);
            try
                for i = 1:3
                    spike_VIP{i} = find(act_VIP(:,i)>0)'.*0.05;
                end
            end

            %SOM
            act_SOM = M(1:t_length,5:7);
            spike_SOM = cell(size(act_SOM,2)+length(spike_VIP),1);
            try
                for i = 1:3
                    spike_SOM{i+length(spike_VIP)} = find(act_SOM(:,i)>0)'.*0.05;
                end
            end
            for i = 1:length(spike_VIP)
                spike_SOM{i} = 80000;
            end

            %PV
            act_PV = [M(1:t_length,8),M(1:t_length,8),M(1:t_length,8)];
            spike_PV = cell(size(act_PV,2)+length(spike_SOM),1);
            try
                for i = 1:3
                    spike_PV{i+length(spike_SOM)} = find(act_PV(:,i)>0)'.*0.05;
                end
            end
            for i = 1:length(spike_SOM)
                spike_PV{i} = 80000;
            end

            %ECS
            %act_ECS = M(1:t_length,9:18);
            act_ECS = M(1:t_length,9);
            spike_ECS = cell(size(act_ECS,2)+length(spike_PV),1);
            try
                for i = 1:size(act_ECS,2)
                    spike_ECS{i+length(spike_PV)} = find(act_ECS(:,i)>0)'.*0.05;
                    if isempty(spike_ECS{i+length(spike_PV)}) == 1
                        spike_ECS{i+length(spike_PV)} = 80000;
                    end
                end
            end
            for i = 1:length(spike_PV)
                spike_ECS{i} = 80000;
            end

            %F
            act_F = M(1:t_length,19);
            spike_F = cell(size(act_F,2)+length(spike_ECS),1);
            try
                for i = 1:size(act_F,2)
                    spike_F{i+length(spike_ECS)} = find(act_F(:,i)>0)'.*0.05;
                    if isempty(spike_F{i+length(spike_ECS)}) == 1
                        spike_F{i+length(spike_ECS)} = 80000;
                    end
                end
            end
            for i = 1:length(spike_ECS)
                spike_F{i} = 80000;
            end

            act_ECS_rest = M(1:t_length,10);
            spike_ECS_rest = cell(size(act_ECS_rest,2)+length(spike_F),1);
            try
                for i = 1:size(act_ECS_rest,2)
                    spike_ECS_rest{i+length(spike_F)} = find(act_ECS_rest(:,i)>0)'.*0.05;
                    if isempty(spike_ECS{i+length(spike_F)}) == 1
                        spike_ECS{i+length(spike_PV)} = 80000;
                    end
                end
            end
            for i = 1:length(spike_F)
                spike_ECS_rest{i} = 80000;
            end

            act_F_rest = M(1:t_length,20);
            spike_F_rest = cell(size(act_F_rest,2)+length(spike_ECS_rest),1);
            try
                for i = 1:size(act_F_rest,2)
                    spike_F_rest{i+length(spike_ECS_rest)} = find(act_F_rest(:,i)>0)'.*0.05;
                    if isempty(spike_F_rest{i+length(spike_ECS_rest)}) == 1
                        spike_F_rest{i++length(spike_ECS_rest)} = 80000;
                    end
                end
            end
            for i = 1:length(spike_ECS_rest)
                spike_F_rest{i} = 80000;
            end


            figure('Renderer', 'painters', 'Position', [10 10 800 500])
            LineFormat = struct();
            LineFormat.Color = b;
            LineFormat.LineWidth = 0.3;
            plotSpikeRaster(spike_VIP,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);
            hold on

            LineFormat.Color = p;
            plotSpikeRaster(spike_SOM,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            LineFormat.Color = g;
            plotSpikeRaster(spike_PV,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            LineFormat.Color = a;
            plotSpikeRaster(spike_ECS,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            LineFormat.Color = o;
            plotSpikeRaster(spike_F,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            LineFormat.Color = a;
            plotSpikeRaster(spike_ECS_rest,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            LineFormat.Color = o;
            plotSpikeRaster(spike_F_rest,'VertSpikePosition',-0.5,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 5000],'LineFormat',LineFormat);

            figure('Renderer', 'painters', 'Position', [10 10 800 200])
            plot(M(:,1),tmp(:,1))
            xlim([0 5000])

            tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
            figure('Renderer', 'painters', 'Position', [10 10 800 200])
            plot(M(:,1),tmp(:,1))
            xlim([0 5000])
            ylim([0.01 0.05])            

        end

        gAMPAECSF(k).v(:,trial+1) = tmp(:,1);
        gAMPAFVIP1(k).v(:,trial+1) = tmp(:,2);
        gAMPAFVIP2(k).v(:,trial+1) = tmp(:,3);
        gAMPAFVIP3(k).v(:,trial+1) = tmp(:,4);

        %M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);
        clear tt

        tt = 0:0.05:80000;
        tt(end) = [];

    end
end


%% Fig. 5D
path = '';% /Volumes/TOSHIBA/fearconditioning/heterogeneous/';
filename = 'het';

trialCOUNT = 40;

for k = 1:length(experiment_)
    experiment = experiment_(k);
    for trial = 0:trialCOUNT-1
        [k trial]

        tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);

        gAMPAECSF(k).v(:,trial+1) = tmp(:,1);%(:,2);%(:,2);
        gAMPAFVIP1(k).v(:,trial+1) = tmp(:,2);
        gAMPAFVIP2(k).v(:,trial+1) = tmp(:,3);
        gAMPAFVIP3(k).v(:,trial+1) = tmp(:,4);

    end
end

M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);
%M(800000:end,:) = [];

tmp_gAMPAECSF(1).v = gAMPAECSF(1).v(1:2000:end,:);
tmp_gAMPAFVIP1(1).v = gAMPAFVIP1(1).v(1:2000:end,:);

figure,
[lineOut(1), fillOut(1)] = stdshade(tmp_gAMPAECSF(1).v',0.2,'r',M(1:2000:end,1)/1000);%,alpha,acolor,F,smth);
xlabel('t [s]')
ylabel('gAMPA ECS-F')
set(gca, 'fontsize', 14)  
ylim([0 0.22])
xlim([0 40])
