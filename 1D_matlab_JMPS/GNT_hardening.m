clearvars;
load GNT1.dat;
load GNT3.dat;
load GNT4.dat;
load GNT5.dat;
purple   = [0.494 0.184 0.556];    % line color purple
blue   = [0 0.447 0.741];        % line color blue
green   = [0.466 0.674 0.188];    % line color green
red   = [1 0 0];              % line color red


noel = 80;         % number of element
dy   = 2.60417e-6*0.2*12/11.2;               % element size along y-direction
cl   = green;


dt   = 0.001;        % integration time step
% 0 2.5 3.5 4.5 7;
tmax = 3;

% material property
E    = 117E3;       % Young's modulus in MPa
m    = 0.001;        % rate sensitivity
e0   = 0.001;       % reference strain rate
e_a  = 0.001;       % appiled strain rate


% define array of stress, strain, grad of strain and Lapacian of strain
stress  = zeros(1,noel);
elastic_strain  = zeros(1,noel);
plastic_strain  = zeros(1,noel);
plastic_dot = zeros(1,noel);
TD_strain = zeros(1,noel);

stress_dot = E*e_a*zeros(1,noel);
grad    = zeros(1,noel);
Laplace = zeros(1,noel);
% sig_Y_local   = 450-225*sin((1:noel)/(noel+1)*pi/2);  % sig_Y in MPa
% sig_Y_local((noel/2):noel) = 225+0*sig_Y_local((noel/2:noel));
sig_Y_local = 450-225*(1:noel)/noel;
sig_Y_tot     = sig_Y_local;
hardening     = zeros(1,noel);

% parameter for gradient hardening
l1   = 1.953e4;           % length parameter for gradient correction
l2   = 0;            % length parameter for Laplace correction
e1   = 0.015;
n1   = 0.6;
e2   = 1E-4;
n2   = 2.0;
h0   = 2E3;

% overall S-S curve
e = 100*(0:dt:tmax)*e_a;
s = 0*e;
h = 0*e;


for i=1:length(e)
% gradient calculation
    for j=2:noel-1
        grad(j) = (plastic_strain(j+1)-plastic_strain(j-1))/2/dy;
    end
    grad(1)=2*grad(2)-grad(3);
    grad(noel)=2*grad(noel-1)-grad(noel-2);
% Laplace calculation
    for j=2:noel-1
        Laplace(j) = (grad(j+1)-grad(j-1))/2/dy;
    end
    Laplace(1) = 2*Laplace(2)-Laplace(3);
    Laplace(noel)=2*Laplace(noel-1)-Laplace(noel-2);
    
    for j=1:noel
        hardening(j) = h0*(1+(sqrt(l1*abs(grad(j))))/(1+(plastic_strain(j)/e2)^n2))/(1+(plastic_strain(j)/e1)^n1);
    end
    

    plastic_dot = e0*(((1-TD_strain).^2).*stress.*(sig_Y_tot).^(-1)).^(1/m);
        
    sig_Y_tot = sig_Y_tot+hardening.*plastic_dot*dt;
    
    plastic_strain = plastic_strain + plastic_dot*dt;
    elastic_strain=(e_a - plastic_dot)*dt;
    stress_dot = E*(e_a - plastic_dot);
    stress = stress + stress_dot*dt;
    TD_strain = -0.3*elastic_strain-0.5*plastic_strain;
    
    s(i) = sum(stress)/noel;
    h(i) = sum(hardening)/noel;
end



hold on

x = (1:noel)*dy;
% 
% figure(1)
% plot(e,s,'Color',cl,'LineWidth',3)
% plot(GNT1(:,1),GNT1(:,2),'--','LineWidth',3);
% plot(GNT3(:,1),GNT3(:,2),'--','LineWidth',3);
% plot(GNT4(:,1),GNT4(:,2),'--','LineWidth',3);
% plot(GNT5(:,1),GNT5(:,2),'--','LineWidth',3);

% title('Stress-Strain curve');
% xlabel('Strain (%)');
% ylabel('Stress (MPa)');
% set(gca,'FontSize',18,'LineWidth',3);
% xlim([0,10]);

% plot(log(1+e/100)*100,h,'Color',cl,'LineWidth',3);
% title('Strain hardening');
% xlabel('True strain (%)');
% ylabel('Work hardening rate (MPa)');
% set(gca,'FontSize',18,'LineWidth',3);
% set(gca,'xtick',0:2:8);
% set(gca,'ytick',300:300:1800);
% xlim([0,8]);
% ylim([300,1800]);

%gradient strengthening
% figure(2)
plot((1:noel)/noel,sig_Y_tot,'Color',cl,'LineWidth',3);
title('Strength profile');
xlabel('Distance');
ylabel('Strength (MPa)');
set(gca,'FontSize',18,'LineWidth',3);
set(gca,'xtick',0:0.2:1);
% % 
% figure(3)
% plot(grad,'LineWidth',3);
% title('Gradient profile');
% xlabel('Distance');
% ylabel('Plastic strain gradient');
% set(gca,'FontSize',18,'LineWidth',3);
% % 
% %figure(4)
% %plot(Laplace,'LineWidth',3);
% %set(gca,'FontSize',18,'LineWidth',3);
% 
% figure(5)
% plot((1:noel)/noel,plastic_strain*100,'LineWidth',3);
% title('Plastic strain profile');
% xlabel('Distance');
% ylabel('Plastic strain (%)');
% set(gca,'FontSize',18,'LineWidth',3);
% set(gca,'xtick',0:0.2:1);
% 
% figure(6)
% plot(TD_strain,'LineWidth',3);
% title('Transverse strain');
% xlabel('Distance');
% ylabel('Strain (MPa)');
% set(gca,'FontSize',18,'LineWidth',3);

% plot((1:noel)/noel,stress_dot/e0/1000,'LineWidth',3);
% plot((1:noel)/noel,hardening/1000,'LineWidth',3);
% xlabel('Distance');
% ylabel('Hardeing (GPa)');
% set(gca,'FontSize',18,'LineWidth',3);
% set(gca,'xtick',0:0.2:1);
% box on
% legend('Tangent modulus','Hardening','Location','northwest');
% ylim([0,900]);
s(length(s))

data = [(1:noel)'/noel/8 sig_Y_tot' plastic_strain'*100];
