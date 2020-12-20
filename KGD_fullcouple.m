% Written By: Wencheng Jin, Idaho National Laboratory (2020)
% Website: https://sites.google.com/view/wenchengjin/software
% Email: wencheng.jin@inl.gov

close all; clear;

PaceParalleltoolbox_r2016b('cores',8)

global CRACK PROP CONTROL

% addpath Full_Couple_solver_constantFracPressure
addpath Full_Coupled_solver

% w = warning('query','last');
% id = w.identifier;
% warning('off',id);


%Mechanical and hydraulic material properties


PROP.E  = 15.96E9;            %Pa
PROP.nu = 0.2;
PROP.eps_cr = 3.5E-5;
PROP.internal_length = 0.05;    %m
PROP.plane_thickness = 1;       %m
PROP.GI = 90;                  %N/m
PROP.GII = PROP.GI;                  %N/m

f_t = PROP.E * PROP.eps_cr;
beta = 1;
PROP.B = 2*PROP.E*beta*PROP.internal_length*f_t/(2*PROP.E*PROP.GI-f_t^2*beta*PROP.internal_length);

PROP.sigmaMax1 = 1e6;              %N/m^2  Pa
PROP.sigmaMax2 = 1e6;
PROP.sigmaMax = 1e6;
PROP.tauMax = 1e6;              %N/m^2  Pa
PROP.lambdaN = 0.05;
PROP.lambdaT = 0.05;
PROP.alpha = 4;
PROP.beta = 4;
PROP.m = PROP.alpha*(PROP.alpha-1)*PROP.lambdaN^2/(1-PROP.alpha*PROP.lambdaN^2);
PROP.n = PROP.beta*(PROP.beta-1)*PROP.lambdaT^2/(1-PROP.beta*PROP.lambdaT^2);
PROP.deltaN = PROP.GI/PROP.sigmaMax*PROP.alpha*PROP.lambdaN*(1-PROP.lambdaN)^(PROP.alpha-1)*(PROP.alpha/PROP.m+1)*(PROP.lambdaN*PROP.alpha/PROP.m+1)^(PROP.m-1);
PROP.deltaT = PROP.GII/PROP.tauMax*PROP.alpha*PROP.lambdaT*(1-PROP.lambdaT)^(PROP.beta-1)*(PROP.beta/PROP.n+1)*(PROP.lambdaT*PROP.beta/PROP.n+1)^(PROP.n-1);
PROP.PenaltyStiffness = 1e8*PROP.sigmaMax/PROP.deltaN;
PROP.dGnt = 0;
PROP.dGtn = 0;
PROP.deltaN_conj = PROP.deltaN-PROP.deltaN*(PROP.dGnt/PROP.GI)^(1/PROP.alpha);
PROP.deltaT_conj = PROP.deltaT-PROP.deltaT*(PROP.dGtn/PROP.GI)^(1/PROP.beta);
PROP.GammaN = -PROP.GI*(PROP.alpha/PROP.m)^PROP.m;
PROP.GammaT = (PROP.beta/PROP.n)^PROP.n ;
PROP.kappa11 = 7e-13;            %m^2   Permeability
PROP.kappa22 = 7e-13;            %m^2
PROP.viscosity = 1.0e-3;        %Pa*s        0.001 Pa*s
PROP.Ks = 36e9;                       %the bulk modulus of skeleton
PROP.BiotAlpha = [0.79 0.79];         %Biot's coefficent, anisotropic due to material properties
PROP.Kf = 3e9;                        %Pa N/m^2  the bulk modulus of water
PROP.porosity = 0.19;                 %20% the porosity of rock(shale)

%Initial cracks: perforated from the boundary of borehole
% CRACK  = [-10 0; -9.9 0];
CRACK  = [0 2.5; 20 2.5];

% read abaqus input file to obtain nodes, coordinates, connectivity, surfaces and sets;
% file = 'KGD.inp';
file = 'multi_stage.inp';

[EXTDISP,BNoset,BElset,Bsurface] = Preprocessor(file);

update=false;  levelSet(update);    % initiate level set value for every node
StateV_initialization(PROP);	      % Initialize problem history dependent varibles
buildNonlocalTable(PROP);       % find gauss points within nonlocal influence zone for each gauss point

% surf={'Surf-2','Surf-3','Surf-4'};
% applied_stress = [0 -1;-2 0;0 1;];
surf={};
applied_stress = [];
[EXTFORCE] = StressBoundary(Bsurface,surf,applied_stress);      % stress distributions

%  set = {'Set-1','Set-2'};
set = {'Set-1'};
pressureB = [0];
[EXTPRESSURE] = PressureBoundary(BNoset,set,pressureB);  %pore presssure boundary

% surf={'Surf-2','Surf-3','Surf-4'};
surf={'Surf-1'};
applied_flux = [0.0];
injection_rate = 0.0002;   % m^2/s
[EXTDISP,EXTFlUX,EXTPRESSURE] = FluxBoundary(Bsurface,surf,applied_flux,injection_rate,EXTDISP,EXTPRESSURE);

CONTROL.Theta = 2/3;
CONTROL.timeI = 0.0;        % Starting time
CONTROL.timeF = 10;        % Ending time  unit second
CONTROL.deltaT = 0.1;         % Time increment
CONTROL.Ncutting = 10;      % Allowed times of cutting back
CONTROL.Niter = 15;         % Maximum number of iteration for each increment
CONTROL.TOL = 1e-5;         % Convergence tolerance
% for simplified model
% NR_solver(CONTROL,PROP,EXTFORCE,EXTDISP,EXTPRESSURE,EXTFlUX,injection_rate)			
NR_solver(CONTROL,PROP,EXTFORCE,EXTDISP,EXTPRESSURE,EXTFlUX)

