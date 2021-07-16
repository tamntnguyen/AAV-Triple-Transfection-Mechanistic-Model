clear; close all; clc
%% Load Data 
cd Data
load E3_totalDNA_mL.mat
dataTotalDNA    = dataTotalDNA_mL;

load E3_mL_cap.mat
dataCap         = dataCap_mL(1:end-2,:);

load E3_gc_B_mL.mat
dataFullCap     = dataGC_mL;

cd ..

% mode = 'presentation';
mode = 'manuscript';

%% Load parameters
load model_params.mat;

%% 
x0      = 7.6e4*[1; 1; 1];              
[t,x]   = simE3(x0,kAll,60,[]);
xCell   = x;

% cell growth
kCD     = [49428.0668639073,0.0112813300989750];
CDfun   = @(t) 1e6 + kCD(1)/kCD(2)*(1 - exp(-kCD(2)*t));
CDvec   = CDfun(t);
figure
plot(t,CDvec);

% Convert to per mL basis
x       = zeros(length(t),size(x,2));
for i = 1:size(x,1)
    x(i,:) = xCell(i,:)*CDvec(i);
end

fracDNA = sum(x(end,[17 19]),2)/sum(x(end,[6 7 8 15 17 19]),2);

%% Plot
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');

ftsz = 25;
line_wdth = 2;
% Replicated DNA vs encapsidated DNA
figure
subplot(1,2,1)
axis square
hold on
plot(t,sum(x(:,[6 7 8 15 17 19]),2),'k','LineWidth',line_wdth);             %Total viral DNA
plot([dataTotalDNA(:,1); dataTotalDNA(:,1); dataTotalDNA(:,1)],...
          [dataTotalDNA(:,2);dataTotalDNA(:,3); dataTotalDNA(:,4)],'k^',...
                        'LineWidth',line_wdth,'MarkerSize',10)
xlabel('Time (hpt)','Interpreter','latex','FontSize',ftsz)

switch mode
    case 'manuscript'
        ylabel('AAV DNA copies/mL','Interpreter','latex','FontSize',ftsz)
    case 'presentation'
        title('AAV DNA copies/mL','Interpreter','latex','FontSize',ftsz)
end

xticks(0:12:60)
set(gca,'FontSize',ftsz)

plot(t,sum(x(:,[17 19]),2),'r','LineWidth',line_wdth);                      % Encapsidated DNA
plot([dataFullCap(:,1); dataFullCap(:,1); dataFullCap(:,1)],...
                [dataFullCap(:,2); dataFullCap(:,3); dataFullCap(:,4)],'rs',...
                            'LineWidth',line_wdth,'MarkerSize',10)
% plot(t,sum(x(:,6:8),2),'m','LineWidth',line_wdth);                          % Plasmid in cell
plot([6 6], [0 7e9], 'k--', 'LineWidth',1.2)                               % Media exchange

ylim([0 7e9])
ylim1=get(gca,'ylim');
xlim1=get(gca,'xlim');

switch mode
    case 'manuscript'
        text(0.1*xlim1(2),0.9*ylim1(2),'\bf{A}','FontSize',26,'Interpreter','latex')
end

% Total capsid vs full capsid
subplot(1,2,2)
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis square
hold on
plot(t,sum(x(:,[16 17 18 19]),2),'b','LineWidth',line_wdth)                 % Total capsid
plot(dataCap(:,1), dataCap(:,2),'bd',...
                'LineWidth',line_wdth,'MarkerSize',10)
        
        
plot(t,sum(x(:,[17 19]),2),'r','LineWidth',line_wdth)                       % Full capsid
plot([dataFullCap(:,1); dataFullCap(:,1); dataFullCap(:,1)],...
                [dataFullCap(:,2); dataFullCap(:,3); dataFullCap(:,4)],'rs',...
                'LineWidth',line_wdth,'MarkerSize',10)
            
% plot(t,sum(x(:,[6 7 8 15 17 19]),2),'k','LineWidth',line_wdth)                                   % Replicated DNA
plot([6 6], [0 6e9], 'k--', 'LineWidth',1.2)
set(gca,'YScale','log','FontSize',ftsz)
ylim([10^8 10^11])

xlim([0 60])
xticks(0:12:60) 

xlabel('Time (hpt)', 'Interpreter','latex','FontSize',ftsz)
switch mode
    case 'manuscript'
        ylabel('Capsids/mL (All Compartments)', 'Interpreter','latex','FontSize',ftsz)
    case 'presentation'
        title('Capsids/mL (All Compartments)', 'Interpreter','latex','FontSize',ftsz)
end



ylim2=get(gca,'ylim');
xlim2=get(gca,'xlim');

switch mode
    case 'manuscript'
        text(0.1*xlim2(2),0.6*ylim2(2),'\bf{B}','FontSize',26,'Interpreter','latex')
end

%%
figure
plot(1,1,'k^',...
     1,1,'k',...)
     1,1,'rs',...
     1,1,'r',...
     1,1,'bd',...
     1,1,'b',...
     'LineWidth',line_wdth,'MarkerSize',10)
leg2 = legend({'Replicated viral DNA - Data',...
                'Replicated viral DNA - Model',...
                'Full virion - Data',...
                'Full virion - Model',...
                'Total capsid - Data',...
                'Total capsid - Model'},'Interpreter','latex',...
                    'FontSize',13,'NumColumns',3);

keyboard

%% Save fig
cd Figures
print -dpsc2 -cmyk Figure4_legend.eps
print -dpsc2 -cmyk Figure4_in-house.eps

axis square