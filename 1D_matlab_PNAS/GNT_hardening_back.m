clearvars;

num_ele = 100;
dy      = 5E-6*1.75/1.75;
e_dot   = 0.001;
dt      = 1e-3;
e0      = 0.001;                    % reference rate
m       = 0.001;                    % strain rate sensitivity

x = (1:num_ele)'/num_ele;
% grain size and twin thickness
d       = 1E-6*(2.5+13.5*x);
l       = 1E-9*(1/30-x*7/360).^(-1);
% for i=1:num_ele/4
%     d(i) = 16E-6;
%     d(i+num_ele/4) = 10E-6;
%     d(i+2*num_ele/4) = 6E-6;
%     d(i+3*num_ele/4) = 2.5E-6;
%     l(i) = 72E-9;
%     l(i+num_ele/4) = 40E-9;
%     l(i+2*num_ele/4) = 38E-9;
%     l(i+3*num_ele/4) = 30E-9;
% end

% material constant
miu = 45E3;
E   = 110E3;
b   = 0.284E-9;
alpha = 0.4;        % Taylor hardening coeff
M   = 2.5;          % Taylor factor

K0  = 2e5;
K1  = 1E10;
K2  = 1.3E3;

C1  = 1e4;
C2  = 1e3;

TIME = 20/dt;
ES  = zeros(TIME,5);

% stress in MPa
stress  = zeros(num_ele,1);
strain  = zeros(num_ele,1);
ep      = zeros(num_ele,1);
grad    = zeros(num_ele,1);
Lap     = zeros(num_ele,1);
back1   = zeros(num_ele,1);
back2   = zeros(num_ele,1);
rho_SSD = 1E13*ones(num_ele,1);      % in SI unit, /m^2
rho_GND = zeros(num_ele,1);
s       = M*alpha*miu*b*rho_SSD.^(0.5);          % slip resistance in MPa
drdy    = zeros(num_ele,1);

for t=1:(20/dt)
    if(t>200/dt)
        e_dot = -0.001;
    end
    if(t>250/dt)
        e_dot = 0.001;
    end
    ES(t,1) = 100*sum(strain)/num_ele;
    ES(t,2) = sum(stress)/num_ele;
    ES(t,3) = sum(s)/num_ele;
    ES(t,4) = sum(back1)/num_ele;
    ES(t,5) = sum(back2)/num_ele;

    eff    = stress-back1-back2;
    ep_dot = e0*sign(eff).*(abs(eff).*(s.^(-1))).^(1/m);
    rho_SSD_dot = (K0*rho_GND + K1*rho_SSD.^(0.5) - K2*rho_SSD).*abs(ep_dot);
    rho_GND_dot = grad/b;
    s_dot  = M*alpha*miu*b*0.5*((rho_SSD).^(-0.5)).*rho_SSD_dot;
    
    strain = ones(num_ele,1)*e_dot*dt + strain;
    stress = stress + E*(e_dot - ep_dot)*dt;
    
    back_sat1 = 0.83*miu*b*l.^(-1);
    back_dot1 = C1*(back_sat1.*ep_dot - back1.*abs(ep_dot));
    
    grad_ave = sum(abs(grad))/num_ele;
    back_sat2 = 11e-5*M*alpha*miu*b*abs(grad/b);
    back_dot2 = C2*(back_sat2 - back2.*sign(ep_dot)).*ep_dot;
    % update ISVs
    s       = s + s_dot*dt;
    back1   = back1 + back_dot1*dt;
    back2   = back2 + back_dot2*dt;
    rho_SSD = rho_SSD + rho_SSD_dot*dt;
%     rho_GND = rho_GND + rho_GND_dot*dt;
    rho_GND = abs(grad/b);
    ep      = ep + ep_dot*dt;
    ep      = smooth(ep);
    for i=2:num_ele-1
        Lap(i) = (ep(i-1)+ep(i+1)-2*ep(i))/dy/dy;
        grad(i) = (ep(i+1)-ep(i-1))/dy/2;
    end
    Lap(1) = 2*Lap(2) - Lap(3);
    Lap(num_ele) = 2*Lap(num_ele-1) - Lap(num_ele-2);
%   Lap(1) = 0;
%   Lap(num_ele) = 0;
    grad(1) = 2*grad(2) - grad(3);
    grad(num_ele) = 2*grad(num_ele-1) - grad(num_ele-2);
%   grad(1) = grad(2);
%   grad(num_ele) = grad(num_ele-1);
    if(t>1)
        h = (ES(t,2)-ES(t-1,2))/abs(e_dot)/dt;
    end
    
    tf = isreal(back1) + isreal(rho_SSD) + isreal(rho_GND);
    if(tf==2)
        disp('Wrong');
        break;
    end
end
h
[ES(t,2) ES(t,4) ES(t,5) ES(t,2)-ES(t,4)-ES(t,5)]
B = [x ep stress rho_SSD rho_GND s back1 back2];
plot(ES(:,1),ES(:,2),ES(:,1),ES(:,3),ES(:,1),ES(:,4),ES(:,1),ES(:,5));
%plot(x,ep);