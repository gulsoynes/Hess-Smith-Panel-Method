clear; clc; close all;

load exp_data.mat       %Loading Experimental data from NACA TR 824 Report
load ThinTheory.mat     %Loading Results from HW4

point = readmatrix('naca1408.txt'); %Loading Airfoil Profile

v = 1.4207e-5;          %Kinematic Viscosity of Air at 10 degree Celcius (m^2/s)
Re = [3e6, 6e6, 9e6];   %Reynold's Numbers
c = 1;                  %Unit Length Chord
V_inf = Re * v / c;     %Velocity with Given Re (m/s)
AoA = -18:1:18;         %Varying Angle of Attack (deg)
N = 200;                %Number of Panel

%Calculating Cl, Cd, Cm with Hess-Smith Panel Method for varying V_inf and
%AoA
for i = 1:length(AoA)
    for j = 1:length(V_inf)
        [Cl(j).a(i,1),Cd(j).a(i,1),Cm(j).a(i,1)] = ...
            hess_panel(point,N,V_inf(j),AoA(i),c);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)   %Plotting Cl versus varying Angle of Attack
plot(cl_Re3(:,1),cl_Re3(:,2),'o-')
hold on
plot(cl_Re6(:,1),cl_Re6(:,2),'s-')
hold on
plot(cl_Re9(:,1),cl_Re9(:,2),'d-')
hold on
plot(AoA, Cl(1).a,'.-')
hold on
plot(AoA, CL_thin)
legend('$Re = 3\times 10^6$ (exp)',...
    '$Re = 6\times 10^6$ (exp)',...
    '$Re = 9\times 10^6$ (exp)',...
    'HSPM',...
    'Thin Airfoil Theory',...
    'Interpreter','Latex')
xlabel('Angle of Attack, $\alpha$, ($deg$)','Interpreter','Latex')
ylabel('Lift Coefficient, $C_L$','Interpreter','Latex')
xlim([-32 32]);
ylim([-2 2]);
xticks([-32:8:32])
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)   %Plotting Cm versus varying Angle of Attack
plot(cm_Re3(:,1),cm_Re3(:,2),'o-')
hold on
plot(cm_Re6(:,1),cm_Re6(:,2),'s-')
hold on
plot(cm_Re9(:,1),cm_Re9(:,2),'d-')
hold on
plot(AoA, Cm(1).a,'.-')
hold on
plot(AoA, cm_ac_thin)
legend('$Re = 3\times 10^6$ (exp)', ...
    '$Re = 6\times 10^6$ (exp)',...
    '$Re = 9\times 10^6$ (exp)', ...
    'HSPM',...
    'Thin Airfoil Theory',...
    'Interpreter','Latex')
xlabel('Angle of Attack, $\alpha$ ($deg$)','Interpreter','Latex')
ylabel('Moment Coefficient, $C_m$','Interpreter','Latex')
xlim([-32 32]);
ylim([-0.5 0.1])
xticks([-32:8:32])
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)   %Plotting Cd versus Cl
plot(cd_cl_Re3(:,1),cd_cl_Re3(:,2),'o-')
hold on
plot(cd_cl_Re6(:,1),cd_cl_Re6(:,2),'s-')
hold on
plot(cd_cl_Re9(:,1),cd_cl_Re9(:,2),'d-')
hold on
plot(Cl(1).a, Cd(1).a,'.-')
legend('$Re = 3\times 10^6$ (exp)', ...
    '$Re = 6\times 10^6$ (exp)', ...
    '$Re = 9\times 10^6$ (exp)', ...
    'HSPM',...
    'Interpreter','Latex')
ylabel('Drag Coefficient, $C_D$','Interpreter','Latex')
xlabel('Lift Coefficient, $C_L$','Interpreter','Latex')
xlim([-1.6 1.6]);
ylim([-0.02 0.02])
xticks([-1.6:0.4:1.6])
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)   %Plotting Cm versus Cl
plot(cm_cl_Re3(:,1),cm_cl_Re3(:,2),'o-')
hold on
plot(cm_cl_Re6(:,1),cm_cl_Re6(:,2),'s-')
hold on
plot(cm_cl_Re9(:,1),cm_cl_Re9(:,2),'d-')
hold on
plot(Cl(1).a, Cm(1).a,'.-')
hold on
plot(CL_thin, cm_ac_thin)
legend('$Re = 3\times 10^6$ (exp)',...
    '$Re = 6\times 10^6$ (exp)', ...
    '$Re = 9\times 10^6$ (exp)', ...
    'HSPM',...
    'Thin Airfoil Theory',...
    'Interpreter','Latex')
ylabel('Moment Coefficient, $C_m$','Interpreter','Latex')
xlabel('Lift Coefficient, $C_L$','Interpreter','Latex')
xlim([-1.6 1.6]);
ylim([-0.5 0.1]);
xticks([-1.6:0.4:1.6])
grid on