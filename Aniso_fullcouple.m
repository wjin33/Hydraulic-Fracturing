close all; clear;

PaceParalleltoolbox_r2016b('cores',6)

global CRACK PROP CONTROL

% addpath Full_Couple_solver_constantFracPressure
addpath Full_Coupled_solver

%Mechanical and hydraulic material properties

PROP.E11  = 20.0E3;            %MPa
PROP.E22  = 10.0E3;            %MPa
PROP.nu12 = 0.2;
PROP.nu23 = 0.2;
PROP.G12 = 15.0E3/2/(1+0.2);   %MPa
PROP.eqeps_1t = 0.9E-4;
PROP.eqeps_2t = 0.8E-4;
PROP.alpha_1t = 4.0E-4;
PROP.alpha_2t = 3.5E-4;
PROP.eqeps_1s = 6.8E-4;
PROP.internal_length = 10;       %mm
PROP.plane_thickness = 1;        %mm
PROP.G1 = 0.19;                  %N/m
PROP.G2 = 0.095;                 %N/m
PROP.sigmaMax1 = 2;              %N/m^2  Pa
PROP.sigmaMax2 = 1;              %N/m^2  Pa
PROP.kappa11 = 4e-10;            %mm^2   Permeability
PROP.kappa22 = 2e-10;            %mm^2
PROP.viscosity = 1.0e-9;         %MPa*s        0.001 Pa*s
PROP.Ks = 36e3;                       %the bulk modulus of skeleton
PROP.BiotAlpha = [0.65 0.75];         %Biot's coefficent, anisotropic due to material properties
PROP.Kf = 3e3;                        %MPa N/mm^2  the bulk modulus of water
PROP.porosity = 0.19;                 %20% the porosity of rock(shale)

%Initial cracks: perforated from the boundary of borehole
% CRACK  = [400 250; 440 250]; % 0 - degree
% CRACK  = [402.6795 240; 437.3205 260]; % 30 - degree
CRACK  = [410 230.6795; 430 269.3205;]; % 60 - degree
% CRACK  = [420 230; 420 270;]; % 90 - degree


% read abaqus input file to obtain nodes, coordinates, connectivity, surfaces and sets;
 file = 'Anisotropy.inp';
% file = 'Anisotropy_stress.inp';

[EXTDISP,BNoset,BElset,Bsurface] = Preprocessor(file);

update=false;  levelSet(update);    % initiate level set value for every node
StateV_initialization(PROP);	      % Initialize problem history dependent varibles
buildNonlocalTable(PROP);       % find gauss points within nonlocal influence zone for each gauss point

% surf={'Surf-2','Surf-3'};
% applied_stress = [0 -4;-2 0;];
surf={};
applied_stress = [];
[EXTFORCE] = StressBoundary(Bsurface,surf,applied_stress);      % stress distributions

set = {'Set-1','Set-2','Set-3','Set-4'};
pressureB = [0 0 0 0];
[EXTPRESSURE] = PressureBoundary(BNoset,set,pressureB);  %pore presssure boundary

% surf={'Surf-2','Surf-3','Surf-4'};
surf={'Surf-1'};
applied_flux = [0.0];
  injection_rate = 20;   % mm^2/s
% injection_rate = 10;   % mm^2/s
[EXTDISP,EXTFlUX,EXTPRESSURE] = FluxBoundary(Bsurface,surf,applied_flux,injection_rate,EXTDISP,EXTPRESSURE);

CONTROL.Theta = 2/3;
CONTROL.timeI = 0.0;        % Starting time
CONTROL.timeF = 0.2;        % Ending time  unit second
CONTROL.deltaT = 0.005;         % Time increment
CONTROL.Ncutting = 10;      % Allowed times of cutting back
CONTROL.Niter = 30;         % Maximum number of iteration for each increment
CONTROL.TOL = 1e-5;         % Convergence tolerance
% for simplified model
% NR_solver(CONTROL,PROP,EXTFORCE,EXTDISP,EXTPRESSURE,EXTFlUX,injection_rate)			
NR_solver(CONTROL,PROP,EXTFORCE,EXTDISP,EXTPRESSURE,EXTFlUX)


