clear all

trialCOUNT = 39; 

path = '';

experiment = 0; %0: full net;
stdp_rule = 1; % 1: asym with larger tau and same amplitude; 

for cases = 1:3
    
    if cases == 1
         external_stimulation = 1; %0: neither CS nor US; 1: only CS; 2: only US; 3: US+CS
         condition = 0; % pre FC
         filename = 'het';

    elseif cases == 2
        external_stimulation = 1;
        condition = 2; % post FC
        filename = 'het';

    elseif cases == 3 % failure
        external_stimulation = 1;
        condition = 2; % post FC
        filename = 'failure_het';

    end

    
    for trial = 0:trialCOUNT-1
        [cases trial]
        
        LFP = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'LFP_',num2str(trial),'.txt']);
        
        M = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'STDP_',num2str(trial),'.txt']);
        
        tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
        
        %{
        if PLOT == 1
        in = 80000;
        fin = 120000;
        figure('Renderer', 'painters', 'Position', [10 10 1100 300])
        plot(LFP(in:fin,1), LFP(in:fin,2),'LineWidth',1.5)
        hold on
        plot(LFP(in:fin,1), LFP(in:fin,5),'LineWidth',1.5)       
        plot(LFP(in:fin,1), LFP(in:fin,4),'LineWidth',1.5)
        plot(LFP(in:fin,1), LFP(in:fin,3),'LineWidth',1.5)
        legend('I AMPA','d-current','p-current','h-current')
        set(gca, 'fontsize', 14)
        end
        %}
        
        % SPECTRAL POWER
        Fs = 20000; % sampling rate
        t = 0:1/Fs:1-1/Fs;
        
        LFP_AMPAcurr = LFP(:,2) + LFP(:,3) + LFP(:,4) + LFP(:,5);
        %LFP(:,2) AMPA;
        %LFP(:,3) h;
        %LFP(:,4) p;
        %LFP(:,5) d;

        LFP_AMPAcurr(1:40000) = []; % remove the first 2000 ms
        
        LFP_VIP = sum(M(:,2:4),2);
        LFP_VIP(1:40000) = []; % remove the first 2000 ms
        
        LFP_PYR = sum(M(:,[11 21]),2);
        LFP_PYR(1:40000) = []; % remove the first 2000 ms
        
        [pxx_AMPAcurr_pmtm(cases).v(:,trial+1),f_AMPAcurr_pmtm(cases).v(:,trial+1)] = pmtm(LFP_AMPAcurr,[],[0:0.1:70],Fs); %pmtm(LFP,[],[0:100],Fs);
        
        %{
        % percentage (like in Davis' paper)
        idx_freq_lowtheta = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)<=6 & f_AMPAcurr_pmtm(cases).v(:,trial+1)>3);
        idx_freq_hightheta = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)<12 & f_AMPAcurr_pmtm(cases).v(:,trial+1)>6);

        lowtheta_power = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_lowtheta,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_lowtheta,trial+1));
        hightheta_power = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_hightheta,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_hightheta,trial+1));

        idx_freq_above3 = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)>=2 & f_AMPAcurr_pmtm(cases).v(:,trial+1)<45);
        total_power = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_above3 ,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_above3 ,trial+1));
        total_power_(cases,trial+1) = total_power;

        %total_power = trapz(f_AMPAcurr_pmtm(cases).v(:,trial+1),pxx_AMPAcurr_pmtm(cases).v(:,trial+1));
        lowtheta_power_perc(cases,trial+1) = lowtheta_power/total_power * 100;
        hightheta_power_perc(cases,trial+1) = hightheta_power/total_power * 100;

        maxlowtheta_power(cases,trial+1) = max(pxx_AMPAcurr_pmtm(cases).v(idx_freq_lowtheta,trial+1));
        maxhightheta_power(cases,trial+1) = max(pxx_AMPAcurr_pmtm(cases).v(idx_freq_hightheta,trial+1));

        %
        idx_freq_between2and12 = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)>=2 & f_AMPAcurr_pmtm(cases).v(:,trial+1)<12);
        total_power_between2and12 = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_between2and12 ,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_between2and12 ,trial+1));

        %total_power = trapz(f_AMPAcurr_pmtm(cases).v(:,trial+1),pxx_AMPAcurr_pmtm(cases).v(:,trial+1));
        lowtheta_power_perc_between2and12(cases,trial+1) = lowtheta_power/total_power_between2and12 * 100;
        hightheta_power_perc_between2and12(cases,trial+1) = hightheta_power/total_power_between2and12 * 100;

        %
        idx_freq_above2 = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)>=2);
        total_power_above2 = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_above2 ,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_above2 ,trial+1));

        %total_power = trapz(f_AMPAcurr_pmtm(cases).v(:,trial+1),pxx_AMPAcurr_pmtm(cases).v(:,trial+1));
        lowtheta_power_perc_above2(cases,trial+1) = lowtheta_power/total_power_above2 * 100;
        hightheta_power_perc_above2(cases,trial+1) = hightheta_power/total_power_above2 * 100;

        lowtheta_power_nonorm(cases,trial+1) = lowtheta_power;
        hightheta_power_nonorm(cases,trial+1) = hightheta_power;

        idx_freq_lowtheta = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)<=6 & f_AMPAcurr_pmtm(cases).v(:,trial+1)>2);
        idx_freq_hightheta = find(f_AMPAcurr_pmtm(cases).v(:,trial+1)<=14 & f_AMPAcurr_pmtm(cases).v(:,trial+1)>6);
        lowtheta_power_wider(cases,trial+1) = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_lowtheta,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_lowtheta,trial+1));
        hightheta_power_wider(cases,trial+1) = trapz(f_AMPAcurr_pmtm(cases).v(idx_freq_hightheta,trial+1),pxx_AMPAcurr_pmtm(cases).v(idx_freq_hightheta,trial+1));
        %}

        
       tmp = readmatrix([path,filename,'_cond',num2str(condition),'_stim',num2str(external_stimulation),'_exp',num2str(experiment),'_rule',num2str(stdp_rule),'AMPA_ECS_PN_',num2str(trial),'.txt']);
        gAMPA(cases).v(trial+1,:,:) = tmp;

        
        %{
        if PLOT == 1
            
            figure('Renderer', 'painters', 'Position', [10 10 1100 300])
            plot(M(:,1),gAMPA(cases).v(trial+1,:,1))
            xlabel('[ms]')
            ylabel('gAMPA')
            set(gca, 'fontsize', 14)
            ylim([0 0.2])

            
        end
        %}
    end
end

f1 = f_AMPAcurr_pmtm(1).v(:,trial+1);
f2 = f_AMPAcurr_pmtm(2).v(:,trial+1);
f3 = f_AMPAcurr_pmtm(3).v(:,trial+1);
pxx1 = mean(pxx_AMPAcurr_pmtm(1).v(:,:),2);
pxx2 = mean(pxx_AMPAcurr_pmtm(2).v(:,:),2);
pxx3 = mean(pxx_AMPAcurr_pmtm(3).v(:,:),2);


figure;
figure('Renderer', 'painters', 'Position', [10 10 800 250])
[lineOut(1), fillOut(1)] = stdshade(pxx_AMPAcurr_pmtm(1).v(:,:)',0.2,'m',f1);
hold on
[lineOut(2), fillOut(2)] = stdshade(pxx_AMPAcurr_pmtm(2).v(:,:)',0.2,'b',f2);
[lineOut(2), fillOut(2)] = stdshade(pxx_AMPAcurr_pmtm(2).v(:,:)',0.2,'k',f2);
xlim([0 45])
ylim([0 4])
legend(lineOut,'pre FC','post FC - learners', 'post FC - non learners')
xlabel('f [Hz]')
ylabel('Spectral Power')
set(gca, 'fontsize', 14)

