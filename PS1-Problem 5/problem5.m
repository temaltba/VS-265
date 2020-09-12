% VS 265 Neural Computation 
% Challenge Problem 5
% Due 09/11/2020
% Tyler Maltba's MATLAB code

clear

%time step
dt = 0.01;

% time domain
t = 0:dt:500;
nt = length(t);

% parameters
vr = -70;    % IC: v(0) = vr = resting potential
C = 100;      % capacitance
Gleak = 5;  % conductivity --> ambient open channels

% Reversal potentials
vna = 55;      % sodium
vk = -92;      % potassium
vcl = -65;     % chloride

%% LIF model with current

% v'(t) = (1/tau)*(vr - v + I(t)/Gleak)
% Runga-Kutta 4

% allocate
v = zeros(1,nt);
v(1) = vr;

% time loop
for i = 2:nt
    
    t0 = t(i-1);

    % No synaptic channels
    % RK4 slopes
    k1 = rhs1(v(i-1), I(t0), vr, C, Gleak);
    k2 = rhs1(v(i-1) + k1*(dt/2), I(t0+dt/2), vr, C, Gleak);
    k3 = rhs1(v(i-1) + k2*(dt/2), I(t0+dt/2), vr, C, Gleak);
    k4 = rhs1(v(i-1) + k3*dt, I(t0+dt), vr, C, Gleak);
    
    % Solution updates
    v(i) = v(i-1) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end

figure(1)
hold on
plot(t,v,'linewidth',2)
set(gca,'fontsize', 22, 'linewidth',2)
xlabel('t ms')
ylabel('V mV')
title('Potential for fixed G = 5 nS')
% legend('C = 50 pF', 'C = 100 pF','C = 200 pF','C = 300 pF')
% title('Potential for fixed C = 100 pF')
% legend('G = 3 nS', 'G = 5 nS','G = 10 nS')

%% Equilibrum solution of Na and K channels

% set v'(t) = +(t) = 0


% Na equilibrium
dGna = 0:.1:50;      % 0 to 50
Gtot = Gleak + dGna;
v = (vr*Gleak + vna*dGna)./Gtot;

% K equilibrium
dGk = dGna;      % 0 to 50
Gtot = Gleak + dGk;
v1 = (vr*Gleak + vk*dGk)./Gtot;

% Na and K equilibrium
Gtot = Gleak + dGna + dGk;
v2 = (vr*Gleak + vna*dGna +vk*dGk)./Gtot;
v22 = (v1+v)/2;

figure(2)
hold on
plot(dGna,v,'-b',dGk,v1,'-r',dGk,v2,'--m',dGk,v22,'--k','linewidth',2)
xline(Gleak,'k')
title('Equilibrium Potentials')
xlabel('\Delta G nS')
ylabel('V mV')
set(gca,'fontsize',20,'linewidth',2)
legend('Na','K','Na and K','Linear','G_{leak}','location','eastoutside')
xlim([0 dGna(end)])

%% Shunting inhibition

dGcl = 10;   
dGna = 0:.1:50;      % from 0 to 50

% Na only
Gtot = Gleak + dGna;
v3 = (vr*Gleak + vna*dGna)./Gtot;

% Cl
Gtot = Gleak + dGcl;
v4 = (vr*Gleak + vcl*dGcl)./Gtot;

% Na and Cl
Gtot = Gleak + dGna + dGcl;
v5 = (vr*Gleak + vna*dGna +vcl*dGcl)./Gtot;
v6 = (v4+v3)/2;

close all
figure(3)
hold on
plot(dGna,v3,'-b',dGna,ones(size(v3))*v4,'-r',dGna,v5,'-m',dGna,v6,'--k','linewidth',2)
xline(Gleak)
set(gca,'fontsize', 20, 'linewidth',2)
xlabel('\Delta G_{Na} nS')
ylabel('V mV')
title('Division vs. Superposition')
legend('Na','Cl','Na and Cl','Linear','G_{leak}','location','eastoutside')
xlim([0 dGna(end)])



function f = rhs1(v1,I1,vr1,C1,Gleak1)
        tau = C1/Gleak1;
        f = (1/tau)*(-v1 + vr1) + I1/C1;
end

function cur = I(tt)

    if tt<100
        cur = 0;
    else
        cur = 100;
    end
end
