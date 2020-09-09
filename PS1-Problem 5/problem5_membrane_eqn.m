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
C = 50;      % capacitance
Gleak = 10;  % conductivity --> ambient open channels

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

% close all
figure(1)
hold on
plot(t,v,'linewidth',2)
set(gca,'fontsize', 22, 'linewidth',2)
xlabel('t [ms]')
ylabel('V [mV]')
ylim([-75 -35])
% title('Potential for fixed G = 5 [nS]')
% legend('C = 50 [pF]', 'C = 100 [pF]','C = 200 [pF]','C = 300 [pF]')
title('Potential for fixed C = 50 [pF]')
legend('G = 3 [nS]', 'C = 5 [nS]','C = 10 [nS]')

%% Equilibrum solution of Na and K channels

% set v'(t) = 0

I1 = zeros(1,nt);
I1(t<100) = 0;
I1(t>=100) = 100;

% Na equilibrium
dGna = 50;      % 0 to 50
Gtot = Gleak + dGna;
v = (vr*Gleak + vna*dGna + I1)/Gtot;

% K equilibrium
dGk = 50;      % 0 to 50
Gtot = Gleak + dGk;
v1 = (vr*Gleak + vk*dGk + I1)/Gtot;

% Na and K equilibrium
% dGk = 50;      % 0 to 50
Gtot = Gleak + dGna + dGk;
v3 = (vr*Gleak + vna*dGna +vk*dGk + I1)/Gtot;

figure(2)
hold on
plot(t,v,'-b',t,v1,'-r',t,v3,'--k','linewidth',2)
% title('Na and K Superpostion')
xlabel('t [ms]')
ylabel('V [mV]')
set(gca,'fontsize',20,'linewidth',2)
legend('\Delta G_{Na} = 50, \Delta G_{k} = 0',...
    '\Delta G_{Na} = 0, \Delta G_{k} = 50',...
    '\Delta G_{Na} = 50, \Delta G_{k} = 50','location','eastoutside')
xlim([0 500])

%% Shunting inhibition

dGcl = 10;   
dGna = 50;      % from 0 to 50
Gtot = Gleak + dGna + dGcl;
tau = C/Gtot;

% allocate
v = zeros(1,nt);
v(1) = vr;

% time loop
for i = 2:nt
    
    t0 = t(i-1);

    % Na and CL synaptic channels
    % RK4 slopes
    k1 = rhs2(v(i-1), I(t0), vr, C, Gleak, vna, vcl, dGna, dGcl, tau);
    k2 = rhs2(v(i-1) + k1*(dt/2), I(t0+dt/2), vr, C, Gleak,...
        vna, vcl, dGna, dGcl, tau);
    k3 = rhs2(v(i-1) + k2*(dt/2), I(t0+dt/2), vr, C, Gleak,...
        vna, vcl, dGna, dGcl, tau);
    k4 = rhs2(v(i-1) + k3*dt, I(t0+dt), vr, C, Gleak,...
        vna, vcl, dGna, dGcl, tau);
    
    % Solution updates
    v(i) = v(i-1) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
end


% close all
figure(3)
hold on
plot(t,v,'--k','linewidth',2)
set(gca,'fontsize', 20, 'linewidth',2)
xlabel('t [ms]')
ylabel('V [mV]')
title('Division vs. Superposition')
% legend('\Delta G_{Na} = 0, \Delta G_{Cl} = 0','\Delta G_{Na} = 0, \Delta G_{Cl} = 10',...
%     '\Delta G_{Na} = 15, \Delta G_{Cl} = 0','\Delta G_{Na} = 15, \Delta G_{Cl} = 10',...
%     '\Delta G_{Na} = 30, \Delta G_{Cl} = 0','\Delta G_{Na} = 30, \Delta G_{Cl} = 10',...
%     '\Delta G_{Na} = 50, \Delta G_{Cl} = 0','\Delta G_{Na} = 50, \Delta G_{Cl} = 10','location','eastoutside')
% legend('\Delta G_{Na} = 0, \Delta G_{Cl} = 10','\Delta G_{Na} = 50, \Delta G_{Cl} = 0',...
%     'Actual: \Delta G_{Na} = 50, \Delta G_{Cl} = 10','SP: \Delta G_{Na} = 15, \Delta G_{Cl} = 10',...
%     'location','eastoutside')


function f = rhs1(v1,I1,vr1,C1,Gleak1)
        tau = C1/Gleak1;
        f = (1/tau)*(-v1 + vr1) + I1/C1;
end

function f2 = rhs2(v1,I1,vr1,C1,Gleak1,vna1,vcl1,dGna1,dGcl1,tau1)

        f2 = -(v1/tau1) + (vr1*Gleak1 +vna1*dGna1 + vcl1*dGcl1 + I1)/C1;
end

function cur = I(tt)

    if tt<100
        cur = 0;
    else
        cur = 100;
    end
end






