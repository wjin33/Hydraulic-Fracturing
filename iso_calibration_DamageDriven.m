close all; clear;
global CRACK PROP CONTROL LocalElement
addpath Mechanical_solver
% PaceParalleltoolbox_r2016b('cores',6)

% addpath Coupled_solver
%Mechanical and hydraulic material properties 

PROP.E  = 15.96E3;            %MPa
PROP.nu = 0.33;

PROP.internal_length = 1;    %mm
PROP.plane_thickness = 1;     %mm
PROP.GI = 0.1;                %N/mm
PROP.GII = PROP.GI;           %N/mm
PROP.sigmaMax = 1E0;             %N/mm^2  MPa
PROP.tauMax = 1E0;               %N/mm^2  MPa
PROP.eps_cr = 1;
% PROP.eps_cr = PROP.sigmaMax/PROP.E;

PROP.lambdaN = 0.01;
PROP.lambdaT = 0.01;
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
PROP.GammaN = -PROP.GI*(PROP.alpha/PROP.m)^PROP.m ;
PROP.GammaT = (PROP.beta/PROP.n)^PROP.n ;
PROP.nonlocal = false;

%Initial cracks: perforated from the boundary of borehole
% CRACK  = [0 0; 323 0 ];
% CRACK  = [170 0; 170 152];
% CRACK  = [];
CRACK  = [0 -1; 0 295];
% read abaqus input file to obtain nodes, coordinates, connectivity, surfaces and sets;
file = 'iso_calibration_coarse.inp';
% file = 'iso_calibration_fine.inp';

[EXTDISP,BNoset,BElset,Bsurface] = Preprocessor(file);

update=false;  levelSet(update);    % initiate level set value for every node
StateV_initialization(PROP);        % Initialize problem history dependent varibles
if PROP.nonlocal
    buildNonlocalTable(PROP,BElset);       % find gauss points within nonlocal influence zone for each gauss point
end

% assign the damage enabled element set
damage_set = 'Set-4';
for i=1:10
    if (strcmp(BElset{i,1},damage_set)>0)
        LocalElement = BElset{i,2};
    end
end


% surf={'Surf-2','Surf-3','Surf-4'};
% applied_stress = [0 -1;-2 0;0 1;];
surf={};
applied_stress = [];
[EXTFORCE] = StressBoundary(Bsurface,surf,applied_stress);      % stress distributions

CONTROL.Theta = 0.6;
CONTROL.timeI = 0.0;        % Starting time
CONTROL.timeF = 1.0;        % Ending time  unit second
CONTROL.deltaT = 0.002;         % Time increment
CONTROL.Ncutting = 15;      % Allowed times of cutting back
CONTROL.Niter =20;         % Maximum number of iteration for each increment
CONTROL.TOL = 1e-5;         % Convergence tolerance
file_name = 'iso_calibration_CZM_';
% file_name = 'iso_calibration_continuum_coarse_';
% file_name = 'iso_calibration_continuum_fine_';
NR_calibration_solver(CONTROL,PROP,EXTFORCE,EXTDISP,file_name)	
