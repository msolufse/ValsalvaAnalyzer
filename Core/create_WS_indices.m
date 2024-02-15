function data = create_WS_indices(patient,plotmarkers)
% function create_WS_indices
% Input: patient name and figure settings
% Output: augmented data strucure and temporary workspace (saved in WS)
% Uses: set_indices to determine VM phases
% Description: generates valsalva maneuver phases 
% create_WS_BP loads BP and raw data file and creates workspace for patient

global Height Sex Weight SPinds DPraw DPinds dt ECGraw ...
       Hraw HRdata PPraw Praw Rdata RRdata SPraw Traw  

%% Tess why do we have Hraw and HRdata (that is confusing)

load(strcat('../WS/',patient,'_WS_temp.mat'))

% Figure properties
fs = plotmarkers.fs;        
lw = round(plotmarkers.lwt/2);
figSize = plotmarkers.figSize;

%% Find Start and Stop 
SPbar = trapz(SPraw(1:1000-1))/(1000 - 1); % was 100 very short to take the average over

f = figure(1); clf; hold on;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient, 'Interpreter', 'none')

subplot(4,1,1); hold on
  plot(Traw,Praw,'Color',[0 0.4470 0.7410],'linewidth',1);
  plot(Traw,SPraw,'b','LineWidth',lw);
  set(gca,'Fontsize',fs);
  title('VM start (min before first peak)','VM end (max before drop & overshoot)');
  ylabel('BP (mmHg)')
  legend('BP','SBP')
  xlim([Traw(1) Traw(end)]);
  ylim([min(Praw)-5 max(Praw)+5])

subplot(4,1,2); hold on
  plot(Traw,HRdata,'b','LineWidth',lw,'Displayname','HR (bpm)');
  set(gca,'Fontsize',fs);
  ylabel('HR (bpm)');
  xlim([0 Traw(end)]);
  ylim([min(HRdata)-2 max(HRdata)+2])

subplot(4,1,3); hold on
  plot(Traw,Rdata,'b','LineWidth',lw);
  set(gca,'Fontsize',fs);
  ylabel('Resp (mV)');
  xlim([0 Traw(end)]);
  
subplot(4,1,4); hold on
  plot(Traw,ECGraw,'b','LineWidth',1);
  set(gca,'Fontsize',fs);
  ylabel('ECG (mV)');
  xlabel('Time (sec)')
  xlim([Traw(1) Traw(end)]);

[val_times, ~ ] = ginput(2);
val_start       = val_times(1);
[~,i_ts_raw]    = min(abs(Traw - val_start));
val_end         = val_times(2);
[~,i_te_raw]    = min(abs(Traw - val_end));

if i_ts_raw >= i_te_raw || isempty(i_te_raw) || isempty(i_ts_raw)
    figure(1); hold on
    title('Error -- click start/stop again')
    [val_times, ~ ] = ginput(2);
    val_start = val_times(1);
    [~,i_ts_raw] = min(abs(Traw - val_start)); % Find the indices of
    % the time data points that are closest to the VM start
    val_end = val_times(2);
    [~,i_te_raw] = min(abs(Traw - val_end)); % Find the indices of
    % the time data points that are closest to the VM end
end

%% Calculate times to cut data before and after VM 

% Calculate available time at start and end of maneuver
timeAvailableS = Traw(i_ts_raw) - Traw(1);
timeAvailableE = Traw(end) - Traw(i_te_raw);

% Rest time before and after
restTimeS = 30;
restTimeE = 30;

% Test if the rest time asked for exceeds available time
if timeAvailableS < restTimeS
    restTimeS = timeAvailableS;
end
if timeAvailableE < restTimeE
    restTimeE = timeAvailableE;
end

% Start and end cut times
timeCutS = timeAvailableS - restTimeS;
timeCutE = timeAvailableE - restTimeE;

% Cut raw data
startTimeV = find(Traw >= Traw(1)   + timeCutS,1);
endTimeV   = find(Traw >= Traw(end) - timeCutE,1);

Traw_cut   = Traw(startTimeV:endTimeV);
Praw_cut   = Praw(startTimeV:endTimeV);
Hraw_cut   = Hraw(startTimeV:endTimeV);
ECGraw_cut = ECGraw(startTimeV:endTimeV);
SPraw_cut  = SPraw(startTimeV:endTimeV);
DPraw_cut  = DPraw(startTimeV:endTimeV);
PPraw_cut  = PPraw(startTimeV:endTimeV);

Rdata_cut  = Rdata(startTimeV:endTimeV);
HRdata_cut = HRdata(startTimeV:endTimeV);
RRdata_cut = RRdata(startTimeV:endTimeV);

% Sample 
v = (1:100:length(Traw_cut))'; 
Tdata    = Traw_cut(v);
ECGdata  = ECGraw_cut(v);
HRdata   = HRdata_cut(v);
Pdata    = Praw_cut(v); 
Rdata    = Rdata_cut(v);
RRdata   = RRdata_cut(v);
SPdata   = interp1(Traw(SPinds),Praw(SPinds),Tdata,'pchip');
DPdata   = interp1(Traw(DPinds),Praw(DPinds),Tdata,'pchip');
PPdata   = SPdata - DPdata;

% Rescale time values to start at 0
val_start = val_start - Tdata(1);
val_end   = val_end - Tdata(1);
Tdata     = Tdata - Tdata(1);
Traw      = Traw - Traw(1);
Traw_cut  = Traw_cut - Traw_cut(1);

i_ts  = find(Tdata >= val_start,1); % index at the beginning of VM
i_t2l = find(Tdata >= val_end,1);   % index at the end of the breath hold

% Set Indices
[phases] = set_indices(Traw_cut,Praw_cut,Tdata,SPdata,HRdata,Rdata,val_start,val_end,patient,plotmarkers);
i_ts       = phases(1); % breath hold start
i_t1       = phases(2); % end of phase I
i_t2e      = phases(3); % end of early phase II
i_t2l      = phases(4); % end of late phase II
i_t3       = phases(5); % end of phase III
i_PRT      = phases(6); % end of pressure recovery time 
i_t4       = phases(7); % end of phase IV

% Baseline Values
i_ts2      = round(i_ts/2);
SPbar      = trapz(SPdata(1:i_ts - 1))/(i_ts - 1);
Hbar       = trapz(HRdata(1:i_ts2))/(i_ts2);
Rbar       = trapz(Rdata(1:i_ts - 1))/(i_ts - 1);
PPbar      = trapz(PPdata(1:i_ts - 1))/(i_ts -1);

% Calculate mean distance between timepoints
dt = mean(diff(Tdata));

% Find max and min HR during respiration
HRdatarespmax = max(HRdata(1:i_ts));
HRdatarespmin = min(HRdata(1:i_ts));
HRdatarespamp = HRdatarespmax - HRdatarespmin;

[HRdatapeaks,~] = findpeaks(HRdata(1:i_ts),'MinPeakProminence',0.25 * HRdatarespamp);
[~,v_locs_data] = findpeaks(-HRdata(1:i_ts),'MinPeakProminence',0.25 * HRdatarespamp);
HRdatavalleys   = HRdata(v_locs_data);

if isempty(HRdatapeaks)
    [HRdatapeaks,~] = max(HRdata(1:i_ts-round(5/dt)));
end
if isempty(HRdatavalleys)
    [HRdatavalleys,~] = min(HRdata(1:i_ts-round(5/dt)));
end

HmaxR = mean(HRdatapeaks);
HminR = mean(HRdatavalleys);

% Make thoracic pressure signal that is strictly positive
Rnew  = 4*(Rdata - mean(Rdata(1:i_ts-1)))/(max(Rdata(1:i_ts-1))-min(Rdata(1:i_ts-1))); % 0 mean value
rhogh = 1060*0.2*9.81/133.322 - 6;
PR    = Rnew + rhogh;

%Thoracic pressure
Pth = zeros(size(Tdata));
for i = 1:length(Tdata)
    t = Tdata(i);
    if (t > Tdata(i_ts)) && (t < Tdata(i_t2l))
        p = 40*(1 - exp(-2*(t - Tdata(i_ts))));
    else
        p = PR(i);
    end
    Pth(i) = p;
end
Pthdata = movmean(Pth,10); 
Pthbar  = trapz(Pthdata(1:i_ts-1))/(i_ts-1); 

%  Find derivatives of pressure data 
    
% Derivative of SBP
dSBPdt = diff(SPdata)/dt;
dSBPdt = [dSBPdt;dSBPdt(end)];

% Derivative of PP
dPPdt = diff(PPdata)/dt;
dPPdt = [dPPdt;dPPdt(end)];

% Mean values before VM
dSBPdt_bar = trapz(dSBPdt(1:i_ts - 1))/(i_ts - 1);
dPPdt_bar  = trapz(dPPdt(1:i_ts - 1))/(i_ts - 1);

% Derivative of Pth
dPthdt = diff(Pthdata)/dt;
dPthdt = [dPthdt;dPthdt(end)];
dPthdt_bar = trapz(dPthdt(1:i_ts - 1))/(i_ts - 1);

% ODE tolerance 
ODE_TOL  = 1e-8;
DIFF_INC = sqrt(ODE_TOL);

gpars.ODE_TOL   = ODE_TOL;
gpars.DIFF_INC  = DIFF_INC;

% Make Data Structure
clear data
data.Traw       = Traw_cut;
data.HRdata     = HRdata_cut;
data.PPraw      = PPraw_cut;
data.Praw       = Praw_cut;
data.ECGraw     = ECGraw_cut;
data.RRdata     = RRdata_cut;

data.Tdata      = Tdata;
data.PPdata     = PPdata;
data.SPdata     = SPdata;
data.HRdata     = HRdata;
data.Pthdata    = Pthdata;
data.Rdata      = Rdata;
data.DPdata     = DPdata;
data.Pdata      = Pdata;
data.RRdata     = RRdata;

data.age        = Age;
data.height     = Height;
data.weight     = Weight;
data.sex        = Sex;

data.patient    = patient;

data.val_start  = val_start;
data.val_end    = val_end;
data.i_ts       = phases(1); % breath hold start
data.i_t1       = phases(2); % end of phase I
data.i_t2e      = phases(3); % end of early phase II
data.i_t2l      = phases(4); % end of late phase II
data.i_t3       = phases(5); % end of phase III
data.i_PRT      = phases(6); % end of pressure recovery time 
data.i_t4       = phases(7); % end of phase IV

data.HminR      = HminR;
data.HmaxR      = HmaxR;
data.Pthbar     = Pthbar;
data.SPbar      = SPbar;
data.PPbar      = PPbar;
data.Hbar       = Hbar;
data.Rbar       = Rbar;

data.dt         = dt;

data.dSBPdt     = dSBPdt;
data.dPPdt      = dPPdt;
data.dPthdt     = dPthdt;
data.dSBPdt_bar = dSBPdt_bar;
data.dPthdt_bar = dPthdt_bar;
data.dPPdt_bar  = dPPdt_bar;

data.gpars = gpars; 

% Four Panel Figure
%X- and Y- axis limits for plots
Tlims   = [Tdata(1) Tdata(end)];
Plims   = [min(Praw)-5 max(Praw)+2];
Hlims   = [min(HRdata)-2 max(HRdata)+2];
Pthlims = [min(Pthdata)-2 max(Pthdata)+2];
ECGlims = 1000*[min(data.ECGraw)*1.1, max(data.ECGraw)*1.1];
dlims   = [-60 60]; 

%Times for VM phases
tsVM  = Tdata(data.i_ts);
t1VM  = Tdata(data.i_t1);
t2eVM = Tdata(data.i_t2e);
t2lVM = Tdata(data.i_t2l);
t3VM  = Tdata(data.i_t3);
t4VM  = Tdata(data.i_t4);

if data.i_t4 <= length(Tdata)
    t4VM = Tdata(data.i_t4);
end

%Identity
I = ones(2,1);

%X values for shaded regions for each VM phase
x1 = [tsVM t1VM t1VM tsVM];
x2 = [t1VM t2lVM t2lVM t1VM];
x3 = [t2lVM t3VM t3VM t2lVM];
x4 = [t3VM t4VM t4VM t3VM];

%Y values for shaded regions for each VM phase
yH   = [Hlims(1) Hlims(1) Hlims(2) Hlims(2)];
yP   = [Plims(1) Plims(1) Plims(2) Plims(2)];
yPth = [Pthlims(1) Pthlims(1) Pthlims(2) Pthlims(2)];
yd   = [dlims(1) dlims(1) dlims(2) dlims(2)];
yECG = [ECGlims(1) ECGlims(1) ECGlims(2) ECGlims(2)];

%Colors
gray  = [.9 .9 .9];
lgray = [.95 .95 .95];

f = figure(3); clf;
set(gcf,'units','points','position',figSize)
f.Units = 'pixels';
sgtitle(patient, 'Interpreter', 'none')

subplot(3,1,1); hold on;
  patch(x1,yP,gray)
  patch(x2,yP,lgray)
  patch(x3,yP,gray)
  patch(x4,yP,lgray)
  plot(t2eVM * I,Plims,'k:','LineWidth',1)
  plot(Tlims, SPbar * I,'k:','LineWidth',1)
  plot(Tlims, PPbar * I,'k:','linewidth',1)
  h1 = plot(data.Traw,data.Praw,'Color',[0 0.4470 0.7410],'LineWidth',1);
  h2 = plot(data.Tdata,data.SPdata,'b','linewidth',lw);
  set(gca,'Fontsize',fs); 
  legend([h1 h2],'BP','SBP','location','southwest')
  ylabel('BP (mmHg)')
  xlim(Tlims);
  ylim(Plims);

subplot(3,1,2); hold on;
  patch(x1,yH,gray)
  patch(x2,yH,lgray)
  patch(x3,yH,gray)
  patch(x4,yH,lgray)
  plot(t2eVM * I,Hlims,'k:','LineWidth',1)
  plot(Tlims, Hbar * I,'k:','LineWidth',1)
  h1 = plot(Tdata, HRdata,'b','LineWidth',lw);
  set(gca,'Fontsize',fs)
  ylabel('HR (bpm)')
  xlim(Tlims)
  ylim(Hlims)

subplot(3,1,3); hold on;
  patch(x1,yPth,gray)
  patch(x2,yPth,lgray)
  patch(x3,yPth,gray)
  patch(x4,yPth,lgray)
  plot(t2eVM * I,Pthlims,'k:','LineWidth',1)
  plot(Tlims, Pthbar*ones(2,1),'k:','linewidth',1)
  plot(Tdata, Pthdata,'Color','b','LineWidth',lw)
  set(gca,'Fontsize',fs); 
  xlabel('Time (s)')
  ylabel('Pth (mmHg)')
  xlim(Tlims)
  ylim(Pthlims)
  

  figFolder = plotmarkers.figFolder;
  fileName = strcat(patient,'_VMphases.png');
  fullFileName = fullfile(figFolder,'WS_data',fileName);
  exportgraphics(f,fullFileName);
   
  bbutton = uicontrol('Parent',f,'Style','pushbutton',...
      'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
      'String','Save and exit','fontsize',fs,'Callback',@closeGenButton);
      
  uiwait(f)

% Save workspace
save(strcat('../WS/',patient,'_WS.mat'),... %Name of file
    'Traw','ECGraw','Hraw','Praw','HRdata','SPraw','DPraw',...
    'RRdata','data','val_start', 'val_end','T_RRint', ...
    'data','SPinds','DPinds', 'patient','Age','Sex','Height','Weight');

end