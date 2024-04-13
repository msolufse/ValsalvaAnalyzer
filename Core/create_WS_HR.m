function [] = create_WS_HR(patient,plotmarkers)
% function create_WS_HR
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses: correct_HR to correct heart rate values
% Description: loads temporary workspace with corrected ECG data calculates
% heart rate and saves temporary workspace with Heart rate output

global raw_data RRint T_RRint sub1 sub2 TR HR_new

% Figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;
if lw > 2
    lw2 = lw-1;
else
    lw2 =1;
end;

%% load temporary workspace and assign vectors
filename = strcat(patient,'_WS_temp.mat');
load(strcat('../WS/',filename),'raw_data','RRint','T_RRint','TR'); 

Traw   = raw_data.Traw;
Traw   = Traw -Traw(1);
ECGraw = raw_data.ECGraw;
Hraw   = raw_data.Hraw;
Praw   = raw_data.Praw;

RRint   = diff(sort(TR));
HR      = 60./RRint;

%% ECG to HR Calculations
f = figure(1); clf; hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
sub1 = subplot(2,1,1); hold on;
    plot(T_RRint,HR,'bo-','linewidth',lw2);
    set(gca,'fontsize',fs)
    ylabel('HR (bpm)')
    xlim([Traw(1),Traw(end)])
    ylim([min(HR)-5 max(HR)+5]);

sub2 = subplot(2,1,2); hold on;
    plot(Traw,ECGraw,'b-')
    set(gca,'fontsize',fs)
    xlabel('Time (s)')
    ylabel('ECG (mV)')
    xlim([Traw(1),Traw(end)])
linkaxes([sub1, sub2],'x')

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_HeartRateECG');
fullFileName = fullfile(figFolder,'WS_data',fileName);
saveas(f,fullFileName,'epsc');

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeButton);

cbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',plotmarkers.posCB,...
    'String','Correct Heart Rate','fontsize',fs,'Callback',@correctButton);

uiwait(f)

% Interpolate over step function and evaluate at Tdata
HRdata = interp1(T_RRint,HR_new,Traw,'pchip');

% Use HR to find RR intervals
RRdata = 60./HRdata;

s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'HRdata','RRdata','-append');

function [] = closeButton(a,~)
    f = a.Parent;
    HR_new = HR;
    close(f)
end % function close botton %

function [] = correctButton(a,~) 
    
    [HR_new,answer_HR] = correct_HR(HR,Traw,T_RRint,patient,plotmarkers);
    
    f = figure(1); clf; hold on;
    set(gcf,'units','points','position',figSize)
    f.Units = 'pixels';
    sgtitle(patient,'Fontsize',fs+2,'FontWeight','bold','Interpreter','none')
    sub1 = subplot(2,1,1); hold on;
        plot(T_RRint,HR_new,'bo-','linewidth',lw2);
        set(gca,'fontsize',fs)
        ylabel('HR (bpm)')
        xlim([Traw(1),Traw(end)])
        ylim([min(HR)-5 max(HR)+5]);
    sub2 = subplot(2,1,2); hold on;
        plot(Traw,ECGraw,'b-')
        set(gca,'fontsize',fs)
        xlabel('Time (s)')
        ylabel('ECG (mV)')
        xlim([Traw(1),Traw(end)])
        
    linkaxes([sub1, sub2],'x')

    figFolder = plotmarkers.figFolder;
    fileName = strcat(patient,'_HeartRateECG');
    fullFileName = fullfile(figFolder,'WS_data',fileName);
    saveas(f,fullFileName,'epsc');

    bbutton = uicontrol('Parent',f,'Style','pushbutton',...
    'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
    'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
   
end % function correct botton %

end % function create_WS_HR %

