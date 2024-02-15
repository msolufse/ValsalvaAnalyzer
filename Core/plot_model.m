function [] = plot_model(pt_id,Outputs,HR,data,pars,opt_flag,plotmarkers)
% function plot_model.m
% Input: patient name, outputs, heart rate, data structure, model
% parameters, flag, and figure properties
% Output: figures saved to the Figure folder
% Description: plots model results for nominal or optimized parameters

% Figure properties
fs = plotmarkers.fs;
lw = round(plotmarkers.lwt);
figSize = plotmarkers.figSize;

% Unpack data structure 
Tdata   = data.Tdata;
SPdata  = data.SPdata;
PPdata  = data.PPdata; 
Hdata   = data.HRdata;
Pthdata = data.Pthdata;
SPbar   = data.SPbar; 
PPbar   = data.PPbar; 
Hbar    = data.Hbar; 
Pthbar  = data.Pthbar; 
i_ts    = data.i_ts;
i_t1    = data.i_t1;
i_t2l   = data.i_t2l;
i_t2e   = data.i_t2e;
i_t3    = data.i_t3;
i_t4    = data.i_t4;
Traw    = data.Traw;
Praw    = data.Praw; 
patient = data.patient;

Tdata = Tdata - Tdata(1); 
Traw  = Traw - Traw(1);

% Model outputs
T_pb = Outputs(:,6);
T_pr = Outputs(:,7);
T_s  = Outputs(:,8); 
T_s(T_s < 0) = 0; 

% Scaling parameters 
H_pb = exp(pars(24)); 
H_pr = exp(pars(25));
H_s  = exp(pars(26)); 

X_pb = T_pb .* H_pb; 
X_pr = T_pr .* H_pr; 
X_s  = T_s  .* H_s;
    
% Set plot limits

% X- and Y- axis limits for plots
Tlims   = [Tdata(1) Tdata(end)]; 
Plims   = [min(Praw)-10 max(Praw)+10]; 
Hlims   = [min(Hdata)-5 max([Hdata' HR'])+5];
Pthlims = [min(Pthdata)-1 max(Pthdata)+1]; 

effmax  = max([T_pb;T_pr;T_s]);
efflims = [-0.1 , max(effmax*1.05,1)]; 

% Times for VM phases
tsVM  = Tdata(i_ts); 
t1VM  = Tdata(i_t1);
t2eVM = Tdata(i_t2e); 
t2lVM = Tdata(i_t2l); 
t3VM  = Tdata(i_t3); 
t4VM  = Tdata(end);

if i_t4 <= length(Tdata)
    t4VM = Tdata(i_t4);
end

% Identity
I = ones(2,1); 

% X values for shaded regions for each VM phase
x1 = [tsVM t1VM t1VM tsVM];
x2 = [t1VM t2lVM t2lVM t1VM]; 
x3 = [t2lVM t3VM t3VM t2lVM];
x4 = [t3VM t4VM t4VM t3VM]; 

% Y values for shaded regions for each VM phase
yH   = [Hlims(1) Hlims(1) Hlims(2) Hlims(2)]; 
yP   = [Plims(1) Plims(1) Plims(2) Plims(2)];
yPth = [Pthlims(1) Pthlims(1) Pthlims(2) Pthlims(2)]; 
yeff = [efflims(1) efflims(1) efflims(2) efflims(2)]; 

% Colors 
gray  = [.9 .9 .9]; 
lgray = [.95 .95 .95];

% Four Panel Figure
f = figure(1); hold on;
set(gcf,'units','points','position',figSize)
if opt_flag == 0
    sgtitle(strcat(pt_id,'_nominal'),'Interpreter', 'none','fontsize',fs)
elseif opt_flag == 1
    sgtitle(strcat(pt_id,'_weighted_optimization'),'Interpreter', 'none','fontsize',fs)
end
subplot(2,2,1); hold on;
    set(gca,'Fontsize',fs)
    patch(x1,yP,gray)
    patch(x2,yP,lgray)
    patch(x3,yP,gray)
    patch(x4,yP,lgray)
    plot(t2eVM * I,Plims,'k:','LineWidth',1)
    plot(Tlims, SPbar * I,'k:','LineWidth',1)
    h1 = plot(Traw,Praw,'Color',[0 0.4470 0.7410],'LineWidth',1);
    h2 = plot(Tdata,SPdata,'b','linewidth',lw);
    legend([h1 h2],'Praw','SBP','location','southwest')
    ylabel('BP (mmHg)')
    xlim(Tlims)
    ylim(Plims)
    
Ms = length(HR)-length(Tdata)+1;
subplot(2,2,2);hold on;
    set(gca,'Fontsize',fs)
    patch(x1,yH,gray)
    patch(x2,yH,lgray)
    patch(x3,yH,gray)
    patch(x4,yH,lgray)
    plot(t2eVM * I,Hlims,'k:','LineWidth',1);
    plot(Tlims, Hbar * I,'k:','LineWidth',1);
    h1 = plot(Tdata, Hdata,'b','LineWidth',lw);
    h2 = plot(Tdata, HR(Ms:end),'Color','m','LineWidth',lw);
    legend([h1 h2],'Data','Model','location','northwest')
    ylabel('HR (bpm)')
    xlim(Tlims)
    ylim(Hlims)

subplot(2,2,3); hold on;
    set(gca,'Fontsize',fs)
    patch(x1,yPth,gray)
    patch(x2,yPth,lgray)
    patch(x3,yPth,gray)
    patch(x4,yPth,lgray)
    plot(t2eVM * I,Pthlims,'k:','LineWidth',1)
    plot(Tlims,Pthbar*ones(2,1),'k:','linewidth',1)
    plot(Tdata, Pthdata,'Color','b','LineWidth',lw)
    xlabel('Time (s)')
    ylabel('Pth (mmHg)')
    xlim(Tlims)
    ylim(Pthlims)

subplot(2,2,4); hold on;
    set(gca,'Fontsize',fs)
    patch(x1,yeff,gray)
    patch(x2,yeff,lgray)
    patch(x3,yeff,gray)
    patch(x4,yeff,lgray)
    plot(t2eVM * I,efflims,'k:','LineWidth',1)
    h1 = plot(Tdata, T_pb(Ms:end),'m','LineWidth',lw);
    %h2 = plot(Tdata, T_pr,'LineWidth',lw); 
    h3 = plot(Tdata, T_s(Ms:end), 'LineWidth',lw,'color',[ 0.30  0 0.5]);
    xlabel('Time (s)')
    ylabel('Baroreflex')
    %legend([h1 h2 h3],'P','R','S','location','northwest')
    legend([h1 h3],'P','S','location','northwest')
    xlim(Tlims)
    ylim([0 1]);

    figFolder = plotmarkers.figFolder;

    if opt_flag == 1
        fileName = strcat(patient,'_optimized.png');
        fullFileName = fullfile(figFolder,'Model_fits',fileName);
    elseif opt_flag == 0
        fileName = strcat(patient,'_nominal.png');
        fullFileName = fullfile(figFolder,'Model_fits',fileName);
    end
    exportgraphics(f,fullFileName);

    bbutton = uicontrol('Parent',f,'Style','pushbutton',...
       'Position',[25,10,figSize(3)*0.25,figSize(4)*0.05],...
       'String','Save and exit','Fontsize',fs,'Callback',@closeGenButton);

    uiwait(f);

end % function plot_model.m %

