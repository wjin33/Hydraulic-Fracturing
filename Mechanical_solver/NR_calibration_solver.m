% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function NR_calibration_solver(CONTROL,PROP,EXTFORCE,EXTDISP,file_name)
%***********************************************************************
% MAIN PROGRAM FOR HYPERELASTIC/ELASTOPLASTIC ANALYSIS
%***********************************************************************

% Set program parameters
% ITRA -- Maximum number of convergence iterations
% TOL  -- The convergence iteration converges when residual < TOL
% ATOL -- if residual > ATOL,then solution deverges besection starts
% NTOL -- The maximum number of bisection times allowed
% TIMS -- [T_start, T_end, T_inc, start load factor, end load factor]
% NOUT -- Print out file ID

%%
global DISPDD DISPTD FORCE GKF PREFTD NODES xyz_frac

%Initialize the field variables and Determine the total number of degree of freedom
[ NEQ] = FieldInitialization();
DISPTD = PREFTD;
DISPDD = sparse(NEQ,1);
% Initialize global stiffness K and residual vector F
GKF = sparse(NEQ,NEQ);
FORCE = sparse(NEQ,1);

ITRA = CONTROL.Niter;
TOL = CONTROL.TOL;   
NTOL = CONTROL.Ncutting;    % NTOL -- The maximum number of bisection times allowed
Theta = CONTROL.Theta;

TIMEI=CONTROL.timeI;				% Starting time
TIMEF=CONTROL.timeF;				% Ending time
deltaT=CONTROL.deltaT;				% Time increment
CUR1=TIMEI;					% Starting load factor
CUR2=TIMEF;					% Ending load factor
DELTA0 = deltaT;					    % Saved time increment
TIME = TIMEI;						% Starting time
TDELTA = TIMEF - TIMEI;				% Time interval for load step
ITOL = 1;						    % Bisection level
TARY=zeros(NTOL,1);					% Time stamps for bisections
%
% Load increment loop
%----------------------------------------------------------------------
ISTEP = 0; FLAG10 = 1;
while(FLAG10 == 1)					% Solution has been converged, start next increment
  FLAG10 = 0; FLAG20 = 1; FLAG30 = 1;
  % FLAG10 incrmental control
  % FLAG20 Bisection control
  % FLAG30 Iteration
  %
  PREFTD = DISPTD;                           % Store converged displacement
  %
  if(ITOL==1)                              % No bisection
    deltaT = DELTA0;
    TARY(ITOL) = TIME + deltaT;
  else                                  % Recover previous bisection
    ITOL = ITOL-1;                      % Reduce the bisection level
    deltaT = TARY(ITOL)-TARY(ITOL+1);		% New time increment
    TARY(ITOL+1) = 0;				    % Empty converged bisection level
%     ISTEP = ISTEP - 1;				% Decrease load increment
  end
  TIME0 = TIME;					% Save the current time
  %
  % Print results Needed
  if (ISTEP>0)
      % Update stresses and history variables
      UPDATE=true; 
      StateV_output(PROP,UPDATE);
      Postprocessor(ISTEP,EXTDISP,file_name); 
  end  % Converged result, need to output
  %
  if (ISTEP  >= 10000)
      Propagation = true;
      while (Propagation)
        Propagation=PropagationLaw();
        if Propagation
          update = true ;
          levelSet(update);
          StateV_update;
          NEQ = 2*max(max(NODES));                         % Number of degrees of freedom for displacement
          TEMP  = zeros(NEQ,1);                              % The displacement as well as the pore pressure is initially zero;
          TEMP(1:size(PREFTD,1),1) = PREFTD ;
          PREFTD = TEMP;
          DISPTD = TEMP;
          heavDOF = nnz(NODES(:,2));                                                % Define the number of Heavi DOF
          if heavDOF > 0, heaviNodes;     end
        end
      end
  end
  
  TIME = TIME + deltaT;				% Increase time
  ISTEP = ISTEP + 1;
  %
  % Check time and control bisection
  while(FLAG20 == 1)				% Bisection loop start
    FLAG20 = 0;
    if ((TIME-TIMEF)>1E-10)			% Time passed the end time
      if ((TIMEF+deltaT-TIME)>1E-10)		% One more at the end time
        deltaT=TIMEF+deltaT-TIME;			% Time increment to the end
        DELTA0=deltaT;				% Saved time increment
        TIME=TIMEF;					% Current time is the end
      else
        ILOAD=ILOAD+1;				% Progress to next load step
        if(ILOAD>NLOAD)				% Finished final load step
          FLAG10 = 0;				% Stop the program
          break;
        else						% Next load step
          TIME=TIME-deltaT;
          deltaT=TIMS(3,ILOAD);
          DELTA0=deltaT;
          TIME = TIME + deltaT;
          TIMEI = TIMS(1,ILOAD);
          TIMEF = TIMS(2,ILOAD);
          TDELTA = TIMEF - TIMEI;
          CUR1 = TIMS(4,ILOAD);
          CUR2 = TIMS(5,ILOAD);
        end
      end
    end
    %
    % Load factor and prescribed displacements
%     FACTOR = CUR1 + (TIME-TIMEI)/TDELTA*(CUR2-CUR1);
     SDISP = deltaT*EXTDISP(:,2)/TDELTA*(CUR2-CUR1);
    
    %
    % Start convergence iteration
    %------------------------------------------------------------------
    ITER = 0;
    DISPDD = zeros(NEQ,1);
    while(FLAG30 == 1)
      FLAG30 = 0;
      ITER = ITER + 1;
      %
      % Initialize global stiffness K and residual vector F
      GKF = sparse(NEQ,NEQ);
      FORCE = sparse(NEQ,1);
      %
      
      % Assmemble the Force(flux) vector from boundary condition
      % Note we used the simply constant value of boundary condition, if
      % the flux or stress changes with time, it should be adjusted
      if size(EXTFORCE,1)>0
        f_index = EXTFORCE(:,1);
        FORCE(f_index,1) = FORCE(f_index,1) + EXTFORCE(:,2);
      end
  
      % Assemble K and F, become very complex, need to know different part 
      if PROP.nonlocal
        updateDomainBeforeNonlocAverage(PROP);       % Calculate the local equivalent strain and store it in STATEV for nonlocal averaging later               
      end
%       updateDomainBeforeNonlocAverage(PROP);       % Calculate the local equivalent strain and store it in STATEV for nonlocal averaging later               
      UPDATE=false; LTAN=true; NLTAN=false;        % Update the residual/internal force vector as well as the stiffness matrix
      stiffnessMatrix(PROP,LTAN,NLTAN,UPDATE);
%       stiffnessMatrix_parallel(PROP,LTAN,NLTAN,NEQ)
      
      %
      % Check convergence
      if( ITER> 1)
        KNOWNDOF =  EXTDISP(:,1);       % find degree of freedom with applied boundary condition
        React_tot = sum(FORCE(KNOWNDOF).^2);            % Total reaction force
        ALLDOF=1:NEQ;
        FREEDOF=setdiff(ALLDOF,KNOWNDOF);
        Residul_tot = sum(FORCE(FREEDOF).^2);           % Total residual force
        RESN=sqrt(Residul_tot/React_tot);
        
        Pressure_error = 0;
        
        OUTPUT(1, ITER, RESN, Pressure_error, TIME, deltaT)
        
        %
        if (RESN<TOL || ( Residul_tot<TOL))
          FLAG10 = 1;
          break;
        end
        %
        if ((ITER>=ITRA))			% Start bisection
          ITOL = ITOL + 1;
          if(ITOL<=NTOL)
            deltaT = 0.5*deltaT;
            TIME = TIME0 + deltaT;
            TARY(ITOL) = TIME;
            DISPTD = PREFTD;
            fprintf(1,'Not converged. Bisecting load increment %3d\n',ITOL);
          else
            fprintf(1,'\t\t *** Max No. of bisection ***\n'); return;
          end
          FLAG20 = 1;
          FLAG30 = 1;
          break;
        end
      end

      % Prescribed displacement BC
      NDISP=size(EXTDISP,1);
      if NDISP~=0
        FIXEDDOF=EXTDISP(:,1);
        GKF(FIXEDDOF,:)=zeros(NDISP,NEQ);
        GKF(FIXEDDOF,FIXEDDOF)=1E15*eye(NDISP);
        %
        FORCE(FIXEDDOF)=0;
        if ITER==1, FORCE(FIXEDDOF) = 1E15*SDISP; end
      end
      
      % Solve the system equation
      if(FLAG20 == 0)
        SOLN = GKF\FORCE;
        DISPDD = DISPDD + SOLN;
        DISPTD = DISPTD + SOLN;
        FLAG30 = 1;
      else
        FLAG30 = 0;
      end
      if(FLAG10 == 1), break; end
    end 							%20 Convergence iteration
  end 								%11 Bisection
end 								%10 Load increment
%
end
function OUTPUT(FLG, ITER, RESN, Pressure_error, TIME, deltaT)
%*************************************************************************
% Print convergence iteration history
%*************************************************************************
%%
  if FLG == 1
    if ITER>2
      fprintf(1,'%27d %14.5e %14.5e\n',ITER,full(RESN),full(Pressure_error));
    else
      fprintf(1,'\n	Time   TimeInc    Iter	 DispResidual    PresResidual\n');
      fprintf(1,'%10.5f %10.3e %5d %14.5e %14.5e\n',TIME,deltaT,ITER,full(RESN),full(Pressure_error));
    end
  end
end

function heaviNodes
% This function assigns nodes enriched with the Heaviside function as
% either above (+1) or below (-1) the crack.

global NODES PSI

for iNode = 1:size(NODES,1)
    if NODES(iNode,2) ~= 0
        NODES(iNode,3) = sign(PSI(iNode));
    end
end
end