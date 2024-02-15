function [] = create_WS_info(patient,plotmarkers)
% function create_WS_info
% Input: patient name and figure settings
% Output: saved temporary workspace (saved in WS)
% Uses:
% Description: Selects data be analyzed and saves temporary workspace with
% patient characteristics and raw ECG, HR, and BP data.

%% Input age, height, and weight
prompt1 = {'Age','Sex (male (m)/female(f))','Height (cm)','Weight (kg)'};
dlgtitle = patient; 
dims = [1 75];

[answer1] = inputdlg(prompt1,dlgtitle,dims);

% Cancel - return to previous screen
if sum(size(answer1) == [0 0]) == 2
    return
end

% Patient Info
Age = answer1{1};
Sex = answer1{2};
Height = answer1{3};
Weight = answer1{4};

%% Load LabChart ECG and blood pressure (BP) data
filename = strcat(patient,'.mat');
load(strcat('../Labchart/',filename)); % must be .mat file

start_time_of_file = 0;

% Indecies for channels for EKG, HR, BP
channel_inds = [1,2,3]; % Channel 1 is ECG, channel 2 HR, and channel 3 BP

% Assemble data
tstart = 0;
sz = size(dataend);
endtime = dataend(1,1);
cols = sz(2);
if cols>1
    for i=2:cols
        endtime=endtime+(dataend(1,i)-dataend(end,i-1));
    end
end

dat = zeros(endtime,4);
if length(tickrate) ~= 1
    if mean(tickrate) == tickrate(1)
        tickrate = tickrate(1);
    else
        error(strcat('Check tickrate - something strange. Index ',num2str(pt)))
    end
end

dummy_time_vec= tstart:1/tickrate:endtime/tickrate-1/tickrate;
val_dat(:,1) = dummy_time_vec+start_time_of_file;
for j = 1:length(channel_inds)
    alldata=data(datastart(channel_inds(j),1):dataend(channel_inds(j),1));
    for i=2:cols
        alldata=[alldata data(datastart(channel_inds(j),i):dataend(channel_inds(j),i))];
    end
    val_dat(:,j+1) = data(datastart(channel_inds(j)):dataend(channel_inds(j)));
end

% Create vector with time (Traw), ECG (ECGraw), heart rate (Hraw), and BP (Praw)
Traw   = val_dat(:,1)-val_dat(1,1);
ECGraw = val_dat(:,2);
Hraw   = val_dat(:,3);
Praw   = val_dat(:,4);

%% Plot data and select time series data for analysis

% figure properties
fs = plotmarkers.fs;
lw = plotmarkers.lwt;
figSize = plotmarkers.figSize;

f = figure(1); clf;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';

sgtitle(patient, 'Interpreter', 'none')
subplot(3,1,1)
    plot(Traw,ECGraw,'b')
    set(gca,'fontsize',fs)
    title(['Select Data (~20s before and after VM)']) 
    ylabel('ECG (mV)')
    xlim([0 max(Traw)])
subplot(3,1,2)
    plot(Traw,Hraw,'b','linewidth',lw)
    set(gca,'fontsize',fs)
    ylabel('HR (bpm)')
    xlim([0 max(Traw)])
    ylim([min(Hraw)-2 max(Hraw)+2])
subplot(3,1,3)
    plot(Traw,Praw,'b')
    set(gca,'fontsize',fs)
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    xlim([0 max(Traw)])
    ylim([min(Praw)-5 max(Praw)+5])

[x,~] = ginput(2);

if x(1)<0 == 1
    istart = 1; % Special case if x(1) is negative
else
    istart = max(find(Traw<x(1)));
end
iend   = max(find(Traw<x(2)));

% Save selected data
TrawO  = Traw;
Traw   = Traw(istart:iend);
ECGraw = ECGraw(istart:iend);
Hraw   = Hraw(istart:iend);
Praw   = Praw(istart:iend);

raw_data.Traw   = Traw;
raw_data.ECGraw = ECGraw;
raw_data.Hraw   = Hraw;
raw_data.Praw   = Praw;

figure(1); hold on
sgtitle(patient, 'Interpreter', 'none')
subplot(3,1,1); hold on
    plot(Traw,ECGraw,'r')
    set(gca,'fontsize',fs)
    title('LabChart data')
    ylabel('ECG (mV)')
    xlim([0 max(TrawO)])
subplot(3,1,2); hold on
    plot(Traw,Hraw,'r','linewidth',lw)
    ylabel('HR (bpm)')
    set(gca,'fontsize',fs)
    xlim([0 max(TrawO)])
    ylim([min(Hraw)-2 max(Hraw)+2])
subplot(3,1,3); hold on
    plot(Traw,Praw,'r')
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    xlim([0,max(TrawO)])
    ylim([min(Praw)-5 max(Praw)+5])

figFolder = plotmarkers.figFolder;
fileName = strcat(patient,'_dataAnalyzed.png');
fullFileName = fullfile(figFolder,'WS_data',fileName);
exportgraphics(f,fullFileName);

bbutton = uicontrol('Parent',f,'Style','pushbutton',...
      'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
 
uiwait(f)
 
% Save temmporary workspace 
s = strcat('../WS/',patient,'_WS_temp.mat');
save(s,'raw_data','Age','Sex','Height','Weight');

end  % function create_WS_info.m %

