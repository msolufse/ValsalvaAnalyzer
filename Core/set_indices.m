function [phases] = set_indices(Traw,Praw,Tdata,SPdata,Hdata,Rdata,val_start,val_end,patient,plotmarkers)
% function set_indices.m
% Input: time, blood pressure, heart rate, valsalva information, patient
% name, and figure properties
% Output: Valsalva maneuver phases
% Uses:
% Description: Determines the valsalva maneuver phases

close all
fs  = plotmarkers.fs;
lw  = plotmarkers.lwt;
if lw > 2
    lw2 = lw-1;
else
    lw2 = 1;
end;
figSize = plotmarkers.figSize;

% Time Indices - breath hold start (ts) & breath hold end (te) !!!
% Calculate mean distance between timepoints
dt = mean(diff(Tdata));

% Rescale times to start at 0
val_start = val_start - Tdata(1);
val_end   = val_end - Tdata(1);
Tdata     = Tdata - Tdata(1);
Traw      = Traw - Traw(1);

% Find the indices of the time data points that are closest to the VM start
% and end times
[~,i_ts]  = min(abs(Tdata - val_start));
[~,i_t2l] = min(abs(Tdata - val_end));

% Set rest before (r1) and after (r2)
i_r1 = 1;
i_r2 = length(Tdata);

% Find steady-state baseline values up to VM start
SPbar = trapz(SPdata(1:i_ts - 1))/(i_ts - 1); % sp

% Find time for end of phase I
% Max 5 seconds after breath hold starts
% make it the peak before the fall (stroke volume decrease response....

[~,k_HRm] = min(Hdata(i_ts:i_ts+round(5/dt))); 
i_HRm = i_ts + k_HRm - 1; 

l_1 = max(i_ts, i_HRm - round(2/dt)); %Find max b/t start of VM and 2s before HR min in phase I 
[~,k_SBPM] = max(SPdata(l_1:i_HRm)); 
i_t1 = l_1 + k_SBPM - 1; 

% Find time for middle of phase II
% Select the point that which the SBP reaches a minimum during phase II
[~,k_t2e] = min(SPdata(i_t1:i_t2l-round(1/dt))); % subtract 1 second to avoid drop off from breath hold
i_t2e = k_t2e + i_t1 - 1;

% Find end of phase III/start of phase II (min 5 seconds after end of BH
[~,k_t3] = min(SPdata(i_t2l:i_t2l+round(5/dt)));
i_t3   = i_t2l + k_t3 - 1; 

% Check if the pressure at the end of phase III is greater than the mean
N = length(SPdata);
if SPdata(i_t3) > SPbar
    k_PRT = find(SPdata(i_t3+1:i_t3+round(5/dt)) <= SPdata(i_t2l),1,'first');
else
    k_PRT = find(SPdata(i_t3+1:end) > SPbar,1,'first');
end
if isempty(k_PRT) == 1
    i_PRT = N;
else
    i_PRT = k_PRT + i_t3;
end

% Find time for end of phase IV  % CHECK!
% Select the point at which the SBP returns to baseline after the overshoot
k_t4 = find(SPdata(i_PRT+1:end) <= SPbar, 1, 'first');
i_t4 = i_PRT+1 + k_t4;
T4    = i_t4*dt;
MinP4 = 10;
MaxP4 = 30/dt;
if isempty(k_t4) == 1
    k_t4 = find(SPdata(i_PRT+1:end) <= SPbar+5, 1, 'first');
    [~,ind] = min(SPdata(i_PRT+k_t4:i_PRT+k_t4+5/dt));
    k_t4 = k_t4 + ind;
    if isempty(k_t4) == 1
        i_t4 = min(N,i_t3+MaxP4);
    else
        i_t4 = i_PRT + k_t4;
    end
else
    while mean(SPdata(i_t4-1:i_t4))-mean(SPdata(i_t4:i_t4+1)) > 0
        i_t4 = i_t4 + 1;
    end
end
T4D = (i_t4-i_t3)*dt;

if T4D < MinP4
    notdone = 0;
    i = round(MinP4/dt)+i_t3;
    while notdone == 0
        if mean(SPdata(i-1:i))-mean(SPdata(i:i+1)) > 0
            i = i+1;
            if i == N
                notdone = 1;
                T4 = i*dt;
            end;
        else
            notdone = 1;
            T4 = i*dt;
            i_t4 = i;
        end;
    end;
end;

%% Plot results
xt(1) = mean(Tdata(i_ts:i_t1));
xt(2) = mean(Tdata(i_t1:i_t2l));
xt(3) = Tdata(i_t2l+2);
xt(4) = mean(Tdata(i_t3:i_t4));

yt_BP(1:4) = max(SPdata)+5;
yt_HR(1:4) = max(Hdata)+5;
str = {'I','II','III','IV'};

f = figure(2);
set(gcf,'units','points','position',figSize);
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')

subplot(3,1,1); hold on
    p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
    p2 = plot(Tdata,SPdata,'b','LineWidth',lw2,'DisplayName','Systolic BP');
    xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
    xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t2e),'--',' ','LineWidth',1)
    xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
    yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
    text(xt,yt_BP,str,'FontSize',fs)
    set(gca,'Fontsize',fs);
    axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
    ylabel('BP (mmHg)');
    legend([p1 p2],{'BP','SBP'});
    xlim([Tdata(i_r1) Tdata(i_r2)]);
    ylim([min(Praw)-5 max(Praw)+5])

subplot(3,1,2); hold on
    plot(Tdata,Hdata,'b','LineWidth',lw2)
    xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
    xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t2e),'--',' ','LineWidth',1)
    xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
    text(xt,yt_HR,str,'FontSize',fs)
    set(gca,'Fontsize',fs);
    xlim([Tdata(i_r1) Tdata(i_r2)]);
    ylim([min(Hdata)-2 max(Hdata)+2]);
    ylabel('HR (bpm)');
    
subplot(3,1,3); hold on;
    xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
    xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t2e),'--',' ','LineWidth',1)
    xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
    xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
    plot(Tdata,Rdata,'b','LineWidth',lw2)
    set(gca,'Fontsize',fs); 
    xlabel('Time (s)')
    ylabel('Respiration')
    xlim([Tdata(i_r1) Tdata(i_r2)]);

% Correct other time indices
answer = questdlg('Accept indices?','Index Correction','Yes','No ','Yes');

% Handle response
while strcmp(answer,'No ')==1
    prompt = {'Choose index to correct.',...
    'Hold control/command to slect multiple.',...
    'Cancel to accept indices.'};

    indices = {'Phase I start','Early Phase II start','Late Phase II start',...
    'Phase III start','Phase IV start','Phase IV end'};

    [indx,~] = listdlg('PromptString',prompt,...
    'SelectionMode','multiple','ListString',indices,...
    'ListSize',[200,300],'Name','Index Correction',...
    'CancelString',{'Cancel'});
    choice = indices(indx);
        
    for i = 1:length(choice)
        switch choice{i}
            case 'Phase I start'
                subplot(3,1,1); hold on
                sgtitle('Click on start of breath hold (Phase I start)');
                title('- Phase I starts at the minimum before the first peak');
                [tsc, ~ ] = ginput(1);
                if isempty(tsc) == 0
                    [~,i_tsc] = min(abs(Tdata - tsc));
                end
                xline(Tdata(i_tsc),'r-','LineWidth',2)
                subplot(3,1,2); hold on
                xline(Tdata(i_tsc),'r-','LineWidth',2)
                subplot(3,1,3); hold on
                xline(Tdata(i_tsc),'r-','LineWidth',2)
                
                ts = tsc;
                i_ts = i_tsc;
                
            case 'Early Phase II start'
                subplot(3,1,1); hold on
                sgtitle('Click on start of early phase II');
                title(' - Phase II starts when SBP peaks after start of breath hold');
                [t1c, ~ ] = ginput(1);
                if isempty(t1c) == 0
                    [~,i_t1c] = min(abs(Tdata - t1c));
                end
                xline(Tdata(i_t1c),'r-','LineWidth',2)
                subplot(3,1,2); hold on
                xline(Tdata(i_t1c),'r-','LineWidth',2)
                subplot(3,1,3); hold on
                xline(Tdata(i_t1c),'r-','LineWidth',2)
                
                t1 = t1c;
                i_t1 = i_t1c;
                
            case 'Late Phase II start'
                subplot(3,1,1); hold on
                sgtitle('Click on start of late phase II');
                title(' - minimum SBP in phase II');
                [t2ec, ~ ] = ginput(1);
                if isempty(t2ec) == 0
                    [~,i_t2ec] = min(abs(Tdata - t2ec));
                end
                xline(Tdata(i_t2ec),'r--','LineWidth',2)
                subplot(3,1,2); hold on 
                xline(Tdata(i_t2ec),'r--','LineWidth',2)
                subplot(3,1,3); hold on
                xline(Tdata(i_t2ec),'r--','LineWidth',2)
                
                t2e = t2ec;
                i_t2e = i_t2ec;
                
            case 'Phase III start'
                subplot(3,1,1); hold on
                sgtitle('Click on end of breath hold (start of Phase III)');
                title(' - Phase III starts at the peak before transient drop in BP (~15 seconds after start of breath hold)');
                [t2lc, ~ ] = ginput(1);
                if isempty(t2lc) == 0
                    [~,i_t2lc] = min(abs(Tdata - t2lc));
                end
                xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                subplot(3,1,2); hold on
                xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                subplot(3,1,3); hold on
                xline(Tdata(i_t2lc),'r-','LineWidth',lw)
                
                i_t2l = i_t2lc;
                
            case 'Phase IV start'
                subplot(3,1,1); hold on
                sgtitle('Click on start of phase IV');
                title(' - Phase IV end after the overshoot return to baseline');
                [t3c, ~ ] = ginput(1);
                if isempty(t3c) == 0
                    [~,i_t3c] = min(abs(Tdata - t3c));
                end
                xline(Tdata(i_t3c),'r-','LineWidth',2)
                subplot(3,1,2); hold on
                xline(Tdata(i_t3c),'r-','LineWidth',2) 
                subplot(3,1,3); hold on
                xline(Tdata(i_t3c),'r-','LineWidth',2) 
                
                t3 = t3c;
                i_t3 = i_t3c;
            case 'Phase IV end'
                subplot(3,1,1); hold on
                sgtitle('Click on start of phase IV');
                title(' - Phase IV starts at the minimum afer phase III');
                [t4c, ~ ] = ginput(1);
                if isempty(t4c) == 0
                    [~,i_t4c] = min(abs(Tdata - t4c));
                end
                xline(Tdata(i_t4c),'r-','LineWidth',2)
                subplot(3,1,2); hold on
                xline(Tdata(i_t4c),'r-','LineWidth',2) 
                subplot(3,1,3); hold on
                xline(Tdata(i_t4c),'r-','LineWidth',2) 
                
                t4 = t4c;
                i_t4 = i_t4c;
        end
    end
    
    xt(1) = mean(Tdata(i_ts:i_t1));
    xt(2) = mean(Tdata(i_t1:i_t2l));
    xt(3) = Tdata(i_t2l+2);
    xt(4) = mean(Tdata(i_t3:i_t4));

    yt_BP(1:4) = max(SPdata)+5;
    yt_HR(1:4) = max(Hdata)+5;
    str = {'I','II','III','IV'};

    f = figure(2);clf
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs,'FontWeight','bold','Interpreter', 'none')
    
    subplot(3,1,1); hold on
        p1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'DisplayName','BP');
        p2 = plot(Tdata,SPdata,'b','LineWidth',2,'DisplayName','Systolic BP');
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        yline(SPbar,'--','Baseline before','linewidth',1,'FontSize',fs)
        text(xt,yt_BP,str,'FontSize',fs)
        axis([Tdata(i_r1) Tdata(i_r2) min(Praw)-5 max(SPdata)+10]);
        ylabel('BP (mmHg)');
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        set(gca,'Fontsize',fs);

    subplot(3,1,2); hold on
        plot(Tdata,Hdata,'b','LineWidth',lw)
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        text(xt,yt_HR,str,'FontSize',fs)
        ylim([min(Hdata)-5 max(Hdata)+10]);
        ylabel('HR (bpm)');
        set(gca,'Fontsize',fs);

    subplot(3,1,3); hold on;
        set(gca,'Fontsize',fs); 
        xline(Tdata(i_ts),'-',' ','LineWidth',1,'LabelVerticalAlignment', 'top','LabelHorizontalAlignment','right','FontSize',fs)
        xline(Tdata(i_t2l),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t1),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t2e),'--',' ','LineWidth',1)
        xline(Tdata(i_t3),'-',' ','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_t4),'-','','LineWidth',1,'FontSize',fs)
        xline(Tdata(i_PRT),'--','','LineWidth',1,'FontSize',fs)
        plot(Tdata,Rdata,'b','LineWidth',lw)
        xlim([Tdata(i_r1) Tdata(i_r2)]);
        xlabel('Time (s)')
        ylabel('Respiration')

    answer = questdlg('Accept indices?', ...
        'Index Correction', ...
        'Yes','No ','Yes');
    msg_not_equal = 'Indices must be in ascending order';
    
    if i_ts < i_t1 == 0 
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Phase I start > Phase II start','Fix',"modal"));
        answer = 'No ';
    elseif i_t1  < i_t2e == 0
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Phase II start > Middle of Phase II','Fix',"modal"));
        answer = 'No ';
    elseif i_t2e < i_t2l == 0 
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Middle of Phase II > Start of Phase III','Fix',"modal"));
        answer = 'No ';
    elseif i_t2l < i_t3  == 0
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Start of Phase III > Start of Phase IV','Fix',"modal"));
        answer = 'No ';
    elseif i_t3  < i_t4  == 0
        msg_not_equal = 'Indices must be in ascending order';
        uiwait(msgbox(msg_not_equal,'Error',"modal"));
        uiwait(msgbox('Start of Phase IV > End of Phase IV','Fix',"modal"));
        answer = 'No ';
    end
end

%% Output
phases = [i_ts; i_t1; i_t2e; i_t2l; i_t3; i_PRT; i_t4];

end % function set_indices.m %
