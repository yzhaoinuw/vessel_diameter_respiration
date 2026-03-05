clear; %close all;
cd('/mnt/nas/202309_Hashmat')
addpath(genpath('Andy')) % for the function Cal_avg_HeartRate_Yue and Andy's script for respiration
CreateDataList

% colors for plotting
col=[0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840;
0 0 1;
1 0 0;
0 1 0;
0 1 1;
1 0 1;
0.25 0.25 0.25;
0.75 0.75 0.75;
0.5 0.5 0.5;
0 0 1;
1 0 0;
0 1 0;
0 1 1;
1 0 1;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25;
0.25 0.25 0.25];

% if var exists
    % ds = var
% else
    % ds = 1:length

% for ds=1:size(data_list,1)
% for ds=ds2plt(55:end)
%for ds=136:size(data_list,1)
for ds=ds
%for ds=ds
    close all;
    cd(data_list{ds,1})
    disp(data_list{ds,1})
    disp(ds)

    clear t_frames inds_plot t_ecgoffset
    
    load('im.mat')
    load('signals.mat');
    
    % use t_frames if it exists
    if exist('t_frames','var')
        t_vect=t_frames';
        disp("using t frame")
    end
    % align the ECG and imaging
    t_ECG = t_ECG - t_ecgoffset;

    %respiration setup
    avg_respiRate = [];
    fs_ECG = 1./mean(diff(t_ECG));

    %Heart rate setup
    avg_HR_Yue_10min = [];

Fc_svm=[.01 0.3]; %small 
Fc_svm2=[.05 .3];
if strcmp(data_list{ds,3},'KX')
    Fc_c=[1.5 10]; % cardiac motion
elseif contains(data_list{ds,3}, 'dex', IgnoreCase=true)
    Fc_c=[3.5 15];
elseif strcmp(data_list{ds,3},'PPA') || strcmp(data_list{ds,3},'PINP')
    Fc_c=[1.5 15];
else
    Fc_c=[6 15];
end
Fc_noise=15;


fs_median=1/median(diff(t_vect));
segtime=10; % length of segments, in minutes
% timepts=[0 10 25:segtime:t_vect(end)/60-segtime];
timepts=[14 25:segtime:t_vect(end)/60-segtime];
if length(timepts)>3
    timepts_plt=[1 2 3];
else
    timepts_plt=2:length(timepts);
end
% seglength=segtime*60*fs_median;
be=-10:0.2:10; % bin edges for svm histograms;

%diam_temp = diam;
% DLC flag
if data_list{ds, 6}
    diam = diam_dlc;
end
% % change back to diam_temp at end

diam_avg=smoothdata(diam,'movmean',10*60*fs_median);
if strcmp(data_list{ds,1},'/mnt/nas/202309_Hashmat/Nov.2023/M24BC-06112023')
    dthresh=[25 50];
    badind=find(diam<min(dthresh) | diam>max(dthresh));
    diam(badind)=diam_avg(badind);
    disp(['warning: you are excluding ' num2str(length(badind)) ' points'])
end

diam_filt_svm = lowpass(diam-diam_avg, Fc_svm(2), fs_median);
diam_filt_svm2 = bandpass(diam-diam_avg, Fc_svm2, fs_median);
 %diam_filt_svm=bandpass(diam-diam_avg,Fc_svm,fs_median,'ImpulseResponse','iir');
%diam_filt_svm=bandpass(diam-diam_avg,Fc_svm,fs_median);
% padlength=1e4;
% diam_filt_svm=bandpass(padarray(diam,[padlength,0],mean(diam)),Fc_svm,fs_mean,'ImpulseResponse','iir');
% diam_filt_svm=diam_filt_svm(padlength+1:end-padlength);

diam_filt_c=bandpass(diam-diam_avg,Fc_c,fs_median);
diam_filt_noise=highpass(diam-diam_avg,Fc_noise,fs_median);

hr_inst=1./(t_ECG(pk_inds_keep(:,2))-t_ECG(pk_inds_keep(:,1))); % this is different than how Andy does it.

fig_results=figure('WindowState','maximized');

if exist("t_frames","var")
    annotation(fig_results,"textbox",'Position',[0.41, 0.9, 0.2, 0.1],'EdgeColor','none',...
        'String',"ECG and Image are perfectly aliged",'FontSize',14,'HorizontalAlignment','center')
else
    if exist("t_frames_error","var")
        annotation(fig_results,"textbox",'Position',[0.41, 0.9, 0.2, 0.1],'EdgeColor','none',...
        'String',["ECG and Image are NOT aliged: off by " num2str(t_frames_error) ' frames'],'FontSize',14,'HorizontalAlignment','center')
    else
         annotation(fig_results,"textbox",'Position',[0.41, 0.9, 0.2, 0.1],'EdgeColor','none',...
        'String',"ECG and Image are NOT aliged",'FontSize',14,'HorizontalAlignment','center')
    end
end

%plot 3
subplot(333), hh1=plot(t_vect/60, diam*umperpix); xlabel('time (min)'), ylabel('diameter (\mum)')
    hold on, hh2=plot(t_vect/60,(diam_filt_svm+diam_avg)*umperpix);
    hold on, hh3=plot(t_vect/60,(diam_avg)*umperpix);
    

% subplot(332), 
%     % plot(t_vect/60, diam_filt_c*umperpix), xlabel('time (min)'), ylabel({'cardiac (8-15 hz) bandpass', 'fitered diameter (\mum)'}),
%     s=scatter(t_vect/60, diam_filt_c*umperpix); xlabel('time (min)'), ylabel({'cardiac (8-15 hz) bandpass', 'fitered diameter (\mum)'}),
%     s.MarkerFaceColor=col(2,:);
%     s.MarkerFaceAlpha=.01;
%     s.MarkerEdgeColor='none';
%     s.SizeData=20;
%     ylim([-2 2])

subplot(335), spectrogram(diam-diam_avg,round(60*fs_median),[],[],fs_median,'yaxis'), title('diameter'), hold on
    ylim(Fc_c)
    xlab=get(gca,'xlabel');
    if strcmp(xlab.String,'Time (hours)')
        xlim([0 t_vect(end)/3600])
    else
        xlim([0 t_vect(end)/60])
    end
    set(gca,'clim',[-30 -10])
subplot(334), spectrogram(diam-diam_avg,round(10*60*fs_median),[],0:.0001:0.5*fs_median,fs_median,'yaxis'), title('diameter'), hold on
    ylim(Fc_svm)
    set(gca,'clim',[0 20])

amp_noise=zeros(size(timepts));
diam_tp=amp_noise; amp_svm=amp_noise; amp_c=amp_noise; hr=amp_noise; w=amp_noise; bp_svm=amp_noise; bpp_svm=amp_noise; histmap=zeros(size(timepts,2),size(be,2)-1);
amp_svm2 = amp_noise; bp_svm2 = amp_noise; svm_tails = amp_noise; %define initial size?
clear w hr

% LOOP of time bin
for i=1:length(timepts)
%     ind=(round(timepts(i)*60*fs_median)+1):round(timepts(i)*60*fs_median)+seglength;
    ind=find(t_vect>timepts(i)*60,1):find(t_vect>(timepts(i)+segtime)*60,1);
    ind_ECG = find(t_ECG>timepts(i)*60,1):find(t_ECG>(timepts(i)+segtime)*60,1);
    ind_HRpeaks = find(t_ECG(pk_inds_keep(:,1))>timepts(i)*60,1):find(t_ECG(pk_inds_keep(:,1))>(timepts(i)+segtime)*60,1);

    %the 99.9th-95th percentile of the bandpass filtered data
    ind=find(t_vect>timepts(i)*60,1):find(t_vect>(timepts(i)+segtime)*60,1);
    svm_tails(i)=prctile(diam_filt_svm(ind),99.9)-prctile(diam_filt_svm(ind),95);

    if min(abs(i-timepts_plt))==0 % is this one of the time points you want to plot?
%         subplot(336), pwelch(diam(ind)-diam_avg(ind),round(5*60*fs_mean),[],0:.0001:fs_mean/2,fs_mean), hold on

% plot 1
subplot(331), pwelch(diam(ind)-diam_avg(ind),[],[],0:.0001:fs_median/2,fs_median), hold on        
        xlim([0 Fc_svm(end)])
        subplot(332), pwelch(diam(ind)-diam_avg(ind),[],[],[],fs_median), hold on
        xlim(Fc_c)        
    end  
    
    [pxx,f]=pwelch((diam(ind)-diam_avg(ind))*umperpix,[],[],[],fs_median);
    if fs_median>2*Fc_c(2) % if the sampling frequency is larger than 2 times the upper cardiac frequency, try to find the heart rate (hr) and w (measure of the heart rate variability)        
        ind_cardiac=find(f>Fc_c(1) & f<Fc_c(2));
        [~,maxind]=max(pxx(ind_cardiac));
        hr(i)=f(ind_cardiac(maxind)); % peak in the pwelch PSD of the diameter signal
        w(i)=iqr(f(pxx(ind_cardiac)>prctile(pxx(ind_cardiac),98)));
%plot 5
subplot(335),
            if strcmp(xlab.String,'Time (hours)')
                plot(timepts(i)/60,hr(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','w'); hold on
            else
                plot(timepts(i),hr(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','w'); hold on
            end
            h_hr=plot(nan,nan,'ow','DisplayName','peak PSD (diam)');
            legend(h_hr)
%         subplot(339), % I swapped out this plot for the band power one
%         b/c it didn't work that well anyway.
%             plot(timepts(i),w(i),'^','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)), hold on
    else
        hr(i)=nan; w(i)=nan;
    end

    diam_tp(i)=mean(diam(ind))*umperpix;
    amp_svm(i)=iqr(diam_filt_svm(ind))*umperpix;
    amp_svm2(i)=iqr(diam_filt_svm2(ind))*umperpix; % new iqr
    
    amp_c(i)=iqr(diam_filt_c(ind))*umperpix;
    amp_noise(i)=iqr(diam_filt_noise(ind))*umperpix;
    bp_svm(i)=trapz(f(f>Fc_svm(1) & f<Fc_svm(2)),pxx(f>Fc_svm(1) & f<Fc_svm(2)));
    bp_svm2(i)=trapz(f(f>Fc_svm2(1) & f<Fc_svm2(2)),pxx(f>Fc_svm2(1) & f<Fc_svm2(2)));%new band power
    bpp_svm(i)=bp_svm(i)/trapz(f(f>Fc_svm(1)),pxx(f>Fc_svm(1)));

    %calculate respiration rate
    resp_trough = -1*resp(ind_ECG);
    avg_Rate = AvgRespiRate(resp_trough,fs_ECG); %call the function
%     [PKS, LOCS, w, p] = findpeaks(resp_trough,"MinPeakDistance",350,"MinPeakProminence",0.35,"MinPeakHeight",0.2);
%     interval = diff(LOCS);
%     interval_sec = interval / fs_ECG;  
%     inst_rate = 60 ./ interval_sec; 
%     avg_Rate = mean(inst_rate);
    avg_respiRate = [avg_respiRate,avg_Rate];
    %disp(avg_Rate);

    %calculate heart rate
    %initialize mask
    peak_mask = true(1,length(detected_r_peaks));
    pk_inds_keep = int32(pk_inds_keep);
    % Find invalid peaks
    bad_peak_indx = find(~ismember(detected_r_peaks, pk_inds_keep(:,1)));
    % Mark invalid peaks as false so this mask can be used to detect
    % continuous segments of valid peaks
    peak_mask(bad_peak_indx) = false; 
    %call the function
    [avg_HR_Yue_10min(i), med_HR(i), var_HR(i)] = Cal_avg_HeartRate_Yue(ind_ECG,detected_r_peaks,fs_ECG,peak_mask);

    % kim's way of calculating heart rate (not as good as Andy's way
    % because it uses a different threshold for good peaks, uses all the
    % good peaks for calculating the heart rate and just filters based on
    % amplitude and time between the peaks with pk_inds_keep
    HR_avg(i)=mean(hr_inst(ind_HRpeaks));
    HR_med(i)=median(hr_inst(ind_HRpeaks));
    HR_var(i)=iqr(hr_inst(ind_HRpeaks));


%plot 3
subplot(333), hold on, plot(timepts(i),diam_tp(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','k'), hold on     

%plot 6
    % plot average respiration rate in 10 min bin
    % plot heart rate in 10 min bin
    
    subplot(336)
    yyaxis left
    plot(timepts(i), avg_HR_Yue_10min(i),'^--'), hold on
    ylabel("avg heart rate (BPM)")
    yyaxis right
    plot(timepts(i), avg_respiRate(i),'o--'), hold on
    ylabel("avg respiration rate (BPM)")


%plot 8
subplot(338), plot(timepts(i),amp_c(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)), hold on
                    diam_filt_c_ind=diam_filt_c(ind);
                    plot(timepts(i),iqr(diam_filt_c_ind(abs(diam_filt_c_ind)<2))*umperpix,'kx'), hold on

%plot 7
subplot(337), plot(timepts(i),amp_svm(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)), hold on

%plot 9
subplot(339), plot(timepts(i),bp_svm(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:)), hold on

    % look at the SVM pdf
    histmap(i,:)=histcounts(diam_filt_svm(ind)*umperpix,be,'Normalization','pdf');
 
end
subplot(331), legend('14-23 min','25-35 min','35-45 min','Location','northeast')
legend([hh1 hh2 hh3],'original','bandpass filtered','moving average')
subplot(336), title('Respirstion and heart rate'), xlabel('Time(min)');
subplot(338), xlabel('time (min)'), ylabel({'cardiac pulsation (\mum)',['IQR [' num2str(Fc_c) '] Hz']} )
subplot(337), xlabel('time (min)'), ylabel({'SVM pulsation (\mum)',['IQR [' num2str(Fc_svm) '] Hz']})
subplot(339), xlabel('time (min)'), ylabel({'SVM Band Power (\mum^2)',['[' num2str(Fc_svm) '] Hz']})

% plot respiration
    load('signals','resp')
% subplot(336),
% spectrogram(resp,60*1000,[],[],1000,'yaxis'), title('resp')
% ylim([0 8])
% set(gca,'clim',[-20 -5])
% title('respiration')

% plot the noise    
subplot(338), hold on, h8=plot(timepts,amp_noise,':','Color',[0.5 0.5 0.5]); %ylabel('amp noise (\mum)')
legend(h8,'noise','Location','best')
plot(timepts,amp_c,'-k','DisplayName','IQR amp')
legend boxoff

% plot the band power percentage
subplot(339), yyaxis right, plot(timepts,bpp_svm,'--^','Color','r'), ylabel('band power %')
    ax=gca; ax.YAxis(2).Color='r';

% change the color of the plots so that it matches the markers
% subplot(331)
%     h=get(gca,'Child'); 
%     egg=length(timepts_plt):-1:1;
%     for j=1:length(timepts_plt)
%         h(egg(j)).Color=col(timepts_plt(j),:);
%     end
% subplot(332)
%     h=get(gca,'Child'); 
%     for j=1:length(timepts_plt)
%         h(egg(j)).Color=col(timepts_plt(j),:);
%     end

% plot the heart rate and variability -- I'm not actually using any of this
% any more
%if fs_median>2*Fc_c(2)
% %     [s,f,t]=spectrogram(diam-diam_avg,round(60*fs_mean),[],[],fs_mean,'yaxis');
%     [s,f,t]=spectrogram(diam-diam_avg,round(10*60*fs_mean),[],[],fs_mean,'yaxis');
%     ind1=find(f>Fc_c(1) & f<Fc_c(2));
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     % opts.Upper = [Inf Inf 1*sqrt(2)];
%     hr=zeros(size(t)); sig=zeros(size(t));
%     clear egg w98 w99 w95
%     for j=1:size(s,2)
%         egg(ind1)=abs(s(ind1,j));
%         [~,hrind]=max(egg);
%         hr2(j)=f(hrind);
%         x=f(ind1);
%         y=abs(s(ind1,j));
%         w98(j)=iqr(x(y>prctile(y,98)));
%         w95(j)=iqr(x(y>prctile(y,95)));
%         w99(j)=iqr(x(y>prctile(y,99)));
% %         [fitresult, gof] = fit(x(y>prctile(y,98)), y(y>prctile(y,98)), fittype( 'gauss1' ));%,opts);
% %         sig(j)=fitresult.c1/sqrt(2);
%     end
% % sig(sig>.99)=nan;
%     
% %     im=imcrop(abs(s),[1 ind1(1) size(s,2) ind1(end)-ind1(1)]);
% %     T=prctile(im(:),95);
% %     imbw=imbinarize(abs(s),T);
% %     imbw(1:ind1(1),:)=0;
% %     imbw(ind1(end):end,:)=0;
% %     imbw=imfill(imbw,'holes');
% % %     imbw=bwareafilt(imbw,1);
% 
%     subplot(332)
%         yyaxis left, plot(t/60,hr2), ylabel('heart rate (hz)')
%         yyaxis right, plot(t/60,w98), ylabel('heart rate IQR (hz)')

%     subplot(339)
%         plot(timepts,w,':'),ylabel('heart rate IQR (hz)') % I commented
%         out this line, which was the last way I was plotting HR
%         variability, in favor of showing the power band.
%end

% label the plot    
[~,name]=fileparts(data_list{ds,1});
subplot(332), title(name)

% set up png save
figs_dir = '/mnt/nas/202309_Hashmat/figs';
type_drug = [upper(data_list{ds, 4}) '-' upper(data_list{ds, 3})];
if ~exist([figs_dir '/' type_drug], 'dir') % checks for type-drug subfolder
    mkdir([figs_dir '/' type_drug]); % makes type-drug subfolder
end

% save the results
save('results','timepts','Fc_noise','Fc_svm','Fc_c','amp_svm','amp_c',...
    'amp_noise','diam_tp','w','hr','diam_filt_svm','diam_filt_c',...
    'diam_filt_noise','diam_avg','fs_median','bp_svm','bpp_svm','avg_respiRate','avg_HR_Yue_10min',...
    'bp_svm2','Fc_svm2','amp_svm2','diam_filt_svm2', 'svm_tails','var_HR','med_HR'); %new variables
% close 

%% phase averaging
if exist('PhaseAvg.mat','file')
    load('PhaseAvg.mat')
    load('signals','t_ECG')
    subplot(338), yyaxis right, plot(timepts,amp,'^--r','DisplayName','phase avg amp'), ylabel('amp from phase avg (\mum)')
    ax=gca;
    ax.YAxis(2).Color='r';
    subplot(335), 
        if strcmp(xlab.String,'Time (hours)')
             h=plot(timepts/60,1./median_dt,'w--','LineWidth',2,'DisplayName','median heart rate (ECG)');
        else
             h=plot(timepts,1./median_dt,'w--','DisplayName','median heart rate from ECG');
        end
        leg=legend([h h_hr]);
        legend boxoff
    
    figPhaseAvg=figure('Position',[2         617        1878         187]);
    subplot(151)
        hr_ecg=1./(t_ECG(pk_inds_keep(:,2))-t_ECG(pk_inds_keep(:,1)));
        h=histogram2(t_ECG(pk_inds_keep(:,1))/3600,hr_ecg,'DisplayStyle','tile','ShowEmptyBins','on','Normalization','pdf');
        h.BinWidth(2)=.2;
        h.BinWidth(1)=100/3600;
        xlabel('Time (hours)')
        ylabel('heart rate from ECG')
        title('PDF')
        hold on, hp=plot(timepts/60,1./median_dt,'w-o','DisplayName','median heart rate (ECG)');
        leg=legend(hp);
        leg.TextColor='white';
        legend box off
     
    for i=1:size(wallspeed,1)
        subplot(152),plot(wave_t,wave_avg_norm(i,:),'Color',col(i,:)); hold on; xlabel('fraction of cardiac cycle'), ylabel('diameter normalized')
        subplot(154),plot(wave_t,wallspeed(i,:),'Color',col(i,:)), hold on, xlabel('fraction of cardiac cycle'), ylabel('wall speed (\mum/s)')
    end

    subplot(153)
        imagesc(timepts,wave_t,wave_avg_norm')
        hold on, plot(timepts,phase_max,'k--o'), legend('maximum')
        ylabel('fraction of cardiac cycle'), xlabel('time (min)')
        cb=colorbar; cb.Label.String='normalized diameter'; 

    subplot(155), 
        imagesc(timepts,wave_t,wallspeed')
        ylabel('fraction of cardiac cycle'), xlabel('time (min)')
        cb=colorbar; cb.Label.String='wallspeed (\mum/s)'; 
        
    savefig(figPhaseAvg,'PhaseAvg_fig','compact')
    print(figPhaseAvg,[figs_dir '/' type_drug '/' name '_PhaseAvg'],'-dpng')
    
end

savefig(fig_results,'results_fig','compact')
print(fig_results,[figs_dir '/' type_drug '/' name '_results'],'-dpng')


%% Large SVM fluctuations
if ~strcmp(data_list{ds,1},'/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M38NO3h-29112023')
%     svm_cutoff_amp=7.8/umperpix;
%     svm_cutoff_amp=3;
      svm_cutoff_amp=0.05*mean(diam);
    figSVM=figure; 
%     [svm_pks,svm_pk_inds]=findpeaks(diam_filt_svm,'MinPeakHeight',svm_cutoff_amp,'MinPeakProminence',2.6/umperpix);
    [svm_pks,svm_pk_inds]=findpeaks(diam_filt_svm,'MinPeakHeight',svm_cutoff_amp,'MinPeakProminence',svm_cutoff_amp*0.5);
    
    % look for the peak in the original diameter signal, not the filtered
    % one
    svm_pk_inds_diam=nan(size(svm_pk_inds));
    pk_neighborhood_sz=round(1*fs_median); % time in seconds * frames/s to get number of frames for neighborhood
    for i=1:length(svm_pk_inds)
        diam_flat=zeros(size(diam));        
        ind_surrounding=max(1,svm_pk_inds(i)-pk_neighborhood_sz):min(svm_pk_inds(i)+pk_neighborhood_sz,length(diam_flat));
        diam_flat(ind_surrounding)=diam(ind_surrounding);%-diam_filt_c(ind_surrounding);
        [~,svm_pk_inds_diam(i)]=max(diam_flat);
    end
    subplot(326),    
        plot(t_vect/60,(diam-diam_avg)*umperpix), hold on
        plot(t_vect/60,diam_filt_svm*umperpix)
        plot(t_vect(svm_pk_inds)/60,svm_pks*umperpix,'go') 
        plot(t_vect(svm_pk_inds_diam)/60,(diam(svm_pk_inds_diam)-diam_avg(svm_pk_inds_diam))*umperpix,'bo') 
        xlabel('time (min)'), ylabel('diameter (\mum)'), legend('original - moving avg','bandpass filtered','svm peaks (bandpass)','svm peaks (original)')

    
    subplot(321), plot(t_vect(svm_pk_inds)/60,svm_pks*umperpix,'--o'), xlabel('time (min)'), ylabel('peak amplitude (\mum)')
    svm_dt=diff(t_vect(svm_pk_inds));
    subplot(322), plot(t_vect(svm_pk_inds(1:end-1))/60,svm_dt,'--o'), xlabel('time (min)'), ylabel('time between peaks (s)')
    title(name)
    
    % [wave_t,wave_avg,wave_err]=bin_and_plot_avg_waveform_selective(t_vect,diam_filt_svm,t_vect,(diam-diam_avg)*umperpix,svm_pk_inds-300,100,1);
    % [wave_t,wave_avg,wave_err]=bin_and_plot_avg_waveform_selective2(t_vect,diam_filt_svm*umperpix,t_vect,(diam-diam_avg)*umperpix,svm_pk_inds-300,100,1,0);
    % [wave_t,wave_avg,wave_err]=bin_and_plot_avg_waveform_selective2(t_vect,diam_filt_svm*umperpix,t_vect,(diam-diam_avg)*umperpix,svm_pk_inds,100,1,0);
%     inds_plot{1}=find(t_vect(svm_pk_inds)<24*60 & t_vect(svm_pk_inds)>300/fs_median);
%     inds_plot{2}=find(t_vect(svm_pk_inds)>25*60 & t_vect(svm_pk_inds)<49*60);
%     inds_plot{3}=find(t_vect(svm_pk_inds)>49*60);
    
    % if baseline
    %   do 1 group
    %   inds_plot{3}=find(t_vect(svm_pk_inds)>10*60);
    % else
    %   do 3 groups

    if (strcmpi(data_list{ds, 3}, 'baseline'))
        inds_plot{1}=find(t_vect(svm_pk_inds)>14*60);
    else
        inds_plot{1}=find(t_vect(svm_pk_inds)<24*60 & t_vect(svm_pk_inds)>14);
        inds_plot{2}=find(t_vect(svm_pk_inds)>=25*60 & t_vect(svm_pk_inds)<35*60);
        inds_plot{3}=find(t_vect(svm_pk_inds)>=35*60);
    end

    for i=1:size(inds_plot,2)
        if length(inds_plot{i})>5 % make sure there are more than 5 waveforms to average over
            [svm_wave_t(:,i),svm_wave_avg(:,i),svm_wave_err(:,i)]=bin_and_plot_avg_waveform_selective2(t_vect,diam_filt_svm*umperpix,t_vect,(diam-diam_avg)*umperpix,svm_pk_inds_diam(inds_plot{i})-300,linspace(0,50,301),0,0);
            close
            svm_wave_speed(:,i)=gradient(svm_wave_avg(:,i),svm_wave_t(:,i)); % microns/s 
            figure(figSVM), 
            subplot(323), plot(svm_wave_t(:,i),svm_wave_avg(:,i)), ylabel('average waveform (\mum)'), xlabel('time (s)'), hold on
            subplot(324), plot(svm_wave_t(:,i),svm_wave_speed(:,i)), xlabel('time (s)'), ylabel('wave speed (\mum/s)'), hold on
            svm_amp_mean(i)=mean(svm_pks(inds_plot{i}))*umperpix;
%             if i==size(inds_plot,2)
%                 svm_dt_mean(i)=mean(svm_dt(inds_plot{i}(1:end-1)));
%             else
%                 svm_dt_mean(i)=mean(svm_dt(inds_plot{i}));
%             end
            svm_dt_mean(i)=mean(diff(t_vect(svm_pk_inds(inds_plot{i})))); % just recalculate the diff, rather than trying to get it from svm_dt

        else
            svm_wave_t(:,i)=nan(300,1); svm_wave_avg(:,i)=nan(300,1); svm_wave_err(:,i)=nan(300,1);
            svm_wave_speed(:,i)=nan(300,1);
            svm_dt_mean(i)=nan;
            svm_amp_mean(i)=nan;
            subplot(323), plot(svm_wave_t(:,i),svm_wave_avg(:,i)), ylabel('average waveform (\mum)'), xlabel('time (s)'), hold on
            subplot(324), plot(svm_wave_t(:,i),svm_wave_speed(:,i)), xlabel('time (s)'), ylabel('wave speed (\mum/s)'), hold on
        end
    end
    subplot(325), 
    hist1=histogram(diam_filt_svm(t_vect<24*60 & t_vect>14*60)*umperpix,be,'Normalization','pdf','DisplayStyle','stairs'); hold on
    hist2=histogram(diam_filt_svm(t_vect>25*60 & t_vect<35*60)*umperpix,be,'Normalization','pdf','DisplayStyle','stairs');
    hist3=histogram(diam_filt_svm(t_vect>35*60)*umperpix,be,'Normalization','pdf','DisplayStyle','stairs'); 
        xlabel('bandpass filtered diameter (\mum)'), ylabel('pdf'), 
        hold on, plot(svm_cutoff_amp*umperpix*[1 1],ylim,'k--')
    subplot(321), hold on, plot([18 30 36+5],svm_amp_mean,'o-','LineWidth',2), plot([24 24],[svm_cutoff_amp*umperpix 20],'k--'), plot([35 35],[svm_cutoff_amp*umperpix 20],'k--')
    subplot(322), hold on, plot([18 30 36+5],svm_dt_mean,'o-','LineWidth',2), plot([24 24],[0 200],'k--'), plot([35 35],[0 200],'k--')
    
    subplot(323), legend('14-24','25-35','>35')
    save('results','-append','svm_pk_inds','svm_pk_inds_diam','svm_wave_t','svm_wave_avg','svm_wave_err','svm_wave_speed','svm_cutoff_amp','inds_plot','svm_dt_mean','svm_amp_mean','be','histmap')
end

figure, imagesc(timepts,0.5*(be(1:end-1)+be(2:end)),histmap')
ylabel('bin centers (\mum)'), xlabel('time (min)'), cb=colorbar; cb.Label.String='SVM PDF';

% close all
end

% for ds=1:size(data_list,1)
% %     load([data_list{ds,1} '/im'],'umperpix')
%     load([data_list{ds,1} '/results'],'diam_filt_svm')
% %     disp(['umperpix = ' num2str(umperpix)])
%     disp(['umperpix = ' num2str(prctile(diam_filt_svm,1))])
% end
% 
% ds2plt=[2:5 7:9 13:size(data_list,1)];
% figure, t=tiledlayout('flow');
% for ds=ds2plt
%     disp(data_list{ds,1})
%     nexttile    
%     load([data_list{ds,1} '/im'],'umperpix','diam')
%     load([data_list{ds,1} '/results'],'diam_filt_svm')
% %     plot(diam_filt_svm*umperpix),
%     plot((diam-smoothdata(diam,'movmean',10*60*55))*umperpix)
%     [~,name]=fileparts(data_list{ds,1}); title(name)
% end
% % xlabel(t,'time (min)'), ylabel(t,'wall speed (\mum/s)')