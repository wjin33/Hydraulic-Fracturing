% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function  StateV_output(PROP,UPDATE)
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

global CONNEC NODES PSI XYZ STATEV DISPDD DISPTD connec_frac xyz_frac statev_frac LocalElement

% global CONNEC NODES PSI XYZ STATEV DISPDD DISPTD FORCE GKF PREFTD connec_frac xyz_frac statev_frac

for iElem= 1:size(CONNEC,1)
    isLocal = ismember(iElem,LocalElement);
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    
    localD  = [N1*2-1 N1*2 N2*2-1 N2*2 N3*2-1 N3*2 N4*2-1 N4*2];    % Traditional index locations for displacement
    iLocD   = 9;                                                             % Next index location
    
    % Traditional element
    X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2); % Nodal x-coordinates
    Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3); % Nodal y-coordinates     
    xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix
    
    if (HEN == 0)                                                           % Unenriched nodes
        
        Ngp = size(STATEV{iElem},2);
        
        for i = 1:Ngp    
            gp = STATEV{iElem}{i}.natural_coodinates;                       % Gauss points
            W = STATEV{iElem}{i}.gauss_weight;                              % Gauss weights               
            xi = gp(1,1); eta = gp(2,1);                                    % Gauss points
            
            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                    % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];

            [Nxy, detJ] = Shape_Function(xi, eta, xyz);                     % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian             
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;

            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,~] = computeNonlocalEquivalentStrain( nonlocal_table );
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            soft_B = STATEV{iElem}{i}.B;

            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal,soft_B);
            
            % Update plastic variables once convergence is reached within
            % one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
                continue;
            end
        end
            
 
    elseif ( HEN > 0 &&  HEN < 4 )                                                            % Enriched element
        
        Ngp = size(STATEV{iElem},2);

        for i = 1:Ngp
            gp = STATEV{iElem}{i}.natural_coodinates;                 % Gauss points
            W = STATEV{iElem}{i}.gauss_weight;                              % Gauss weights               
            xi = gp(1,1); eta = gp(2,1);                                % Gauss points
            
            [Nxy, detJ] = Shape_Function(xi, eta, xyz);                 % Determinant of the Jacobian and the derivative of shape functions with respect to x/y                                 

            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                     % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];

            Benru = [];
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            index = 1;
            for iN = 1:4
                if NN(iN,2) ~= 0
                    psi1 = PSI(N1);                                         % Psi level set value at node 1
                    psi2 = PSI(N2);                                         % Psi level set value at node 2
                    psi3 = PSI(N3);                                         % Psi level set value at node 3
                    psi4 = PSI(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
    
                    Hgp = sign(psi);                                        % Heaviside value at current gauss point
                    Hi  = NN(iN,3);                                         % Nodal Heaviside value
                    H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value
                    
                    Benru(:,(2*index-1):(2*index)) = [Nxy(1,iN)*H    0;
                                                      0              Nxy(2,iN)*H;
                                                      Nxy(2,iN)*H    Nxy(1,iN)*H];

                    index = index+1;
                    if (i == 1)
                        localD(iLocD:(iLocD+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
                        iLocD = iLocD+2;

                    end
                end    
            end

            Bu = [Bu Benru];
            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;
            
            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,~] = computeNonlocalEquivalentStrain( nonlocal_table );
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            soft_B = STATEV{iElem}{i}.B;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal,soft_B);
                        
            % Update plastic variables once convergence is reached within one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
               continue;
            end
        end

    elseif  ( HEN == 4 )                                                    % Fully enriched element
                                                      
        if ( numel(PSI) == 0 )
            PN = [0 0 0 0]; 
        else
            PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
        end
        [~,~,J] = subDomain(3,PN,xyz);                                  % Full Heaviside enrichment

        Ngp = size(STATEV{iElem},2);

        for i = 1:Ngp
            gp = STATEV{iElem}{i}.natural_coodinates;                    % Gauss points
            W = STATEV{iElem}{i}.gauss_weight;                           % Gauss weights               
            xi = gp(1,1); eta = gp(2,1);                                 % Gauss points
            
            Ji   = [J(i,1) J(i,2);J(i,3) J(i,4)];                        % Jacobian of subdomain
            detJ = det(Ji);                                              % Determinant of the Jacobian
            [Nxy, ~] = Shape_Function(xi, eta, xyz);
            N  = 1/4*[(1-xi)*(1-eta)  (1+xi)*(1-eta)  ...                % Shape functions
                      (1+xi)*(1+eta)  (1-xi)*(1+eta)];
            
            Benru = [];
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            index = 1;
            for iN = 1:4
                if NN(iN,2) ~= 0
                    psi1 = PSI(N1);                                         % Psi level set value at node 1
                    psi2 = PSI(N2);                                         % Psi level set value at node 2
                    psi3 = PSI(N3);                                         % Psi level set value at node 3
                    psi4 = PSI(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
    
                    Hgp = sign(psi);                                        % Heaviside value at current gauss point
                    Hi  = NN(iN,3);                                         % Nodal Heaviside value
                    H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value
                    
                    Benru(:,(2*index-1):(2*index)) = [Nxy(1,iN)*H    0;
                                                     0              Nxy(2,iN)*H;
                                                     Nxy(2,iN)*H    Nxy(1,iN)*H];
                    index = index+1;
                    if (i == 1)
                        localD(iLocD:(iLocD+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
                        iLocD = iLocD+2;
                    end
                end    
            end

            Bu = [Bu Benru];
            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;
            
            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,~] = computeNonlocalEquivalentStrain( nonlocal_table );
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            soft_B = STATEV{iElem}{i}.B;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal,soft_B);
                        
            % Update plastic variables once convergence is reached within one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
                continue;
            end

        end
        
        I_frac = find( connec_frac(:,1) == iElem );
        
        N1_frac  = connec_frac(I_frac,2);                                              % Node 1 for current fracture segment
        N2_frac  = connec_frac(I_frac,3);                                              % Node 2 for current fracture segment
        
%         X1_frac = xyz_frac(N1_frac,4);
%         X2_frac = xyz_frac(N2_frac,4);
        
        tan_frac = xyz_frac(N2_frac,2:3) - xyz_frac(N1_frac,2:3);
        tan_frac = tan_frac./sqrt(sum(tan_frac.^2));
        norm_frac = [-tan_frac(2) tan_frac(1)];
        
        lambda = [tan_frac; norm_frac];
        
        localD_enr  = localD(1,9:16);
        
        Ngp = size(statev_frac{I_frac},2);
        
        for ig = 1:Ngp

            gp_ele = statev_frac{I_frac}{ig}.natural_coordinate_ele;            % Gauss points
            histo_jump = statev_frac{I_frac}{ig}.jump; 
            
            xi = gp_ele(1,1); eta = gp_ele(2,1);                                % Gauss points
            
            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                     % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];
                   

            NT = [N(1)   0     N(2)  0     N(3)  0     N(4)  0 ;
                  0      N(1)  0     N(2)  0     N(3)  0     N(4);];
     
            jump = lambda*NT*DISPTD(localD_enr);
            fracwidth = jump(2);
            
            new_jump=zeros(2,1);
            if jump(1) >= histo_jump(1)
                new_jump(1) = jump(1);
            end
            if jump(2) >= histo_jump(2)
                new_jump(2) = jump(2);
            end

            if UPDATE
                %%%%%%
                statev_frac{I_frac}{ig}.width = fracwidth;
                statev_frac{I_frac}{ig}.jump = new_jump;
                continue;
            end

        end
    end 
end

end

function [stressN, DamageN, kappaN]=Gauss_sig_update(PROP,EquivStrain,strainN,kappa,Damage,isLocal,soft_B)
    % Inputs:
    % PROP = [E0, nu0, epsilon_0, epsilon_f];
    %*D = 4*4 elastic stiffness matrix
    % stress = [s11, s22, s12,   s33];
    % strain = [e11, e22, 2*e12, e33];
    % Plane strain problem e33=0
    %%
    
        if isLocal
            if ( EquivStrain > kappa )
                kappa = EquivStrain;
                Damage = 1-PROP.eps_cr/kappa*exp(-soft_B*(kappa - PROP.eps_cr));
            end
        end
        
        DamageN = Damage;
        kappaN = kappa;
        
        MATC=zeros(4,4);
    
        D= (1-DamageN)*PROP.E/(1+PROP.nu)/(1-2*PROP.nu);
        MATC(1,1)=(1-PROP.nu)*D;
        MATC(1,2)=PROP.nu*D;
        MATC(1,3)=PROP.nu*D;
        MATC(1,4)=0;
        MATC(2,1)=PROP.nu*D;
        MATC(2,2)=(1-PROP.nu)*D;
        MATC(2,3)=PROP.nu*D;
        MATC(2,4)=0;
        MATC(3,1)=PROP.nu*D;
        MATC(3,2)=PROP.nu*D;
        MATC(3,3)=(1-PROP.nu)*D;
        MATC(3,4)=0;
        MATC(4,1)=0;
        MATC(4,2)=0;
        MATC(4,3)=0;
        MATC(4,4)=(1-2*PROP.nu)/2*D;
    
        stressN = MATC*[strainN(1:3,1); 2*strainN(4,1)];    %updated stress
    end

        function [nonlocal_equ_eps,scale] = computeNonlocalEquivalentStrain( nonlocal_table )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            global  STATEV
        
            nonlocal_equ_eps = 0;
            scale = 0;
        
            for i = 1:size(nonlocal_table,1)
                iElem = nonlocal_table(i,1);
                iGP = nonlocal_table(i,2);
                weight = nonlocal_table(i,3);
                volume = STATEV{iElem}{iGP}.volume;
                equivalentEPS = STATEV{iElem}{iGP}.EquivStrain; 
                nonlocal_equ_eps = nonlocal_equ_eps + equivalentEPS.*(weight*volume);
                scale = scale + (weight*volume);
            end
            nonlocal_equ_eps = nonlocal_equ_eps./scale;
            if nonlocal_equ_eps < 0
                stop=1;
            end
        
        end

function [Nxy, detJ] = Shape_Function(xi, eta, Elxy)
%******************************************************************************
% Compute shape function, derivatives, and determinant of 4 Node plane element
%******************************************************************************
%%

 Nxi  = 1/4*[-(1-eta)   1-eta  1+eta  -(1+eta)];          % Derivative of shape functions with respect to x
 Neta = 1/4*[-(1-xi)  -(1+xi)  1+xi       1-xi];          % Derivative of shape functions with respect to y
                
 Jacobi = [Nxi;Neta;]*Elxy;
 
 detJ = det(Jacobi);

 InvJacobi = Jacobi\eye(2);
 
 Nxy = InvJacobi*[Nxi;Neta;];
 
end

function [traction] =  cohesiveLaw(jump, histo_jump)

    global PROP
	delt_max = histo_jump(1);
	deln_max = histo_jump(2);
	delt= abs(jump(1));
	deln= jump(2);
    sign_dt=1;
    if (jump(1) < 0) 
        sign_dt=-1;
    end
    traction = zeros(2,1);
    if ( deln < 0 )
        traction(2)=PROP.PenaltyStiffness*deln;
    elseif( deln >= 0 && deln <= PROP.deltaN  && delt <= PROP.deltaT_conj )
        if (deln >= deln_max )
            traction(2)=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-deln/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m)*...
                                (PROP.GammaT*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n+PROP.dGtn);
        else
            traction(2)=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-deln_max/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln_max/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^PROP.m)*...
                                (PROP.GammaT*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n+PROP.dGtn)*deln/deln_max;
            
        end
    end
    if ( deln >= 0 && deln <= PROP.deltaN_conj && delt <= PROP.deltaT )
        if ( delt >= delt_max )
            traction(1)=(PROP.GammaT/PROP.deltaT)*( (PROP.n*(1-delt/PROP.deltaT)^PROP.beta)*(PROP.n/PROP.beta+delt/PROP.deltaT)^(PROP.n-1)-...
                                PROP.beta*(1-delt/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n)*...
                                (PROP.GammaN*(1-deln/PROP.deltaN)^PROP.alpha*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m+PROP.dGtn)*sign_dt;  
        else
            traction(1)=(PROP.GammaT/PROP.deltaT)*( (PROP.n*(1-delt_max/PROP.deltaT)^PROP.beta)*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^(PROP.n-1)-...
                                PROP.beta*(1-delt_max/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^PROP.n)*...
                                (PROP.GammaN*(1-deln/PROP.deltaN)^PROP.alpha*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m+PROP.dGtn)*delt*sign_dt/delt_max;
        end
    end

end

function [stiffness] =  cohesiveStiffness(jump, histo_jump)

    global PROP
    
	delt_max = histo_jump(1);
	deln_max = histo_jump(2);
	delt= abs(jump(1));
	deln= jump(2);
    sign_dt=1;
    if (jump(1) < 0) 
        sign_dt=-1;
    end
    stiffness = zeros(2,2);
    if ( deln < 0 )
        stiffness(2,2)=PROP.PenaltyStiffness;
        stiffness(2,1)=0;
    elseif( deln >= 0 && deln <= PROP.deltaN  && delt <= PROP.deltaN_conj )
        if (deln >= deln_max )
            stiffness(2,2)=(PROP.GammaN/PROP.deltaN^2)*( ((PROP.m^2-PROP.m)*(1-deln/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^(PROP.m-2)+...
                                ((PROP.alpha^2-PROP.alpha)*(1-deln/PROP.deltaN)^(PROP.alpha-2)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m)-...
                                2*PROP.alpha*PROP.m*(1-deln/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^(PROP.m-1))*...
                                (PROP.GammaT*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n+PROP.dGtn);
            stiffness(2,1)=(PROP.GammaN*PROP.GammaT/PROP.deltaN/PROP.deltaT)*( ( PROP.m*(1-deln/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m)*...
                                ( PROP.n*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^(PROP.n-1)-...
                                  PROP.beta*(1-delt/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n )*sign_dt;
        else
            stiffness(2,2)=(PROP.GammaN/PROP.deltaN)*( (PROP.m*(1-deln_max/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln_max/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^PROP.m)*...
                                (PROP.GammaT*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n+PROP.dGtn)/deln_max;
            stiffness(2,1)=(PROP.GammaN*PROP.GammaT/PROP.deltaN/PROP.deltaT)*( ( PROP.m*(1-deln_max/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln_max/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln_max/PROP.deltaN)^PROP.m)*...
                                ( PROP.n*(1-delt/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt/PROP.deltaT)^(PROP.n-1)-...
                                  PROP.beta*(1-delt/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n )*sign_dt*deln/deln_max;
        end
    end
    if ( deln >= 0 && deln <= PROP.deltaN_conj && delt <= PROP.deltaT )
        if ( delt >= delt_max )
            stiffness(1,1)=(PROP.GammaT/PROP.deltaT^2)*( ((PROP.n^2-PROP.n)*(1-delt/PROP.deltaT)^PROP.beta)*(PROP.n/PROP.beta+delt/PROP.deltaT)^(PROP.n-2)+...
                                ((PROP.beta^2-PROP.beta)*(1-delt/PROP.deltaT)^(PROP.beta-2)*(PROP.n/PROP.beta+delt/PROP.deltaT)^PROP.n)-...
                                2*PROP.beta*PROP.n*(1-delt/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt/PROP.deltaT)^(PROP.n-1))*...
                                (PROP.GammaN*(1-deln/PROP.deltaN)^PROP.alpha*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.n+PROP.dGnt);
            stiffness(1,2)=stiffness(2,1);
        else
            stiffness(1,1)=(PROP.GammaT/PROP.deltaT)*( (PROP.n*(1-delt_max/PROP.deltaT)^PROP.beta)*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^(PROP.n-1)-...
                                PROP.beta*(1-delt_max/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^PROP.n)*...
                                (PROP.GammaN*(1-deln/PROP.deltaN)^PROP.alpha*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m+PROP.dGtn)/delt_max;
            stiffness(1,2) = (PROP.GammaN*PROP.GammaT/PROP.deltaN/PROP.deltaT)*( ( PROP.m*(1-deln/PROP.deltaN)^PROP.alpha)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^(PROP.m-1)-...
                                PROP.alpha*(1-deln/PROP.deltaN)^(PROP.alpha-1)*(PROP.m/PROP.alpha+deln/PROP.deltaN)^PROP.m)*...
                                ( PROP.n*(1-delt_max/PROP.deltaT)^PROP.beta*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^(PROP.n-1)-...
                                  PROP.beta*(1-delt_max/PROP.deltaT)^(PROP.beta-1)*(PROP.n/PROP.beta+delt_max/PROP.deltaT)^PROP.n )*sign_dt*delt/delt_max;    
        end
    end

end

function [gp,gw,J] = subDomain(npt,psi,xyz)
% This function subdivides enriched elements and determines the guass 
% points and weights to be used in the integration during the assembly of 
% the stiffness matrix.
%%
corner = [1 2 3 4 1];
node   = [-1 -1;1 -1;1 1;-1 1];

% Loop through crack level set function
if isempty(psi) == 0
    for i = 1:4
        n1 = corner(i);
        n2 = corner(i+1);
        if psi(n1)*psi(n2) < 0
            r    = psi(n1)/(psi(n1)-psi(n2));
            pnt  = (1-r)*node(n1,:)+r*node(n2,:);
            xi   = pnt(1); eta = pnt(2);
            N    = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                        (1+xi)*(1+eta);(1-xi)*(1+eta)];
            xpnt = dot(N,xyz(1:4,1)');
            ypnt = dot(N,xyz(1:4,2)');
            xyz  = [xyz;xpnt ypnt];
            node = [node;pnt];
        end
    end
end

% Find the triangles defining the subdomains
warning off MATLAB:delaunayn:DuplicateDataPoints
tri = delaunay(node(:,1),node(:,2));

% Loop over subtriangles to get quadrature points and weights
[q,w] = gauss(npt,'TRI');

pt = 1;
for e = 1:size(tri,1)
    coord = node(tri(e,:),:);
    xyzl  = xyz(tri(e,:),:);
    for i = 1:length(w)
        xi = q(i,1); eta = q(i,2);
        N  = [1-xi-eta;xi;eta];
        gp(pt,:) = N'*coord;
        gw(pt,1) = w(i)/2;
        J(pt,:)  = [-xyzl(1,1)+xyzl(2,1) -xyzl(1,2)+xyzl(2,2)...
                    -xyzl(1,1)+xyzl(3,1) -xyzl(1,2)+xyzl(3,2)];
        pt = pt+1;
    end
end

end

function [EquivStrain] = computeEquivalentStrain( strain)
    %%
        [~,strain_principal]= eigs([strain(1,1) strain(4,1); strain(4,1) strain(2,1)]);
        EquivStrain = sqrt(((strain_principal(1,1)+abs(strain_principal(1,1)))/2)^2+...
                           ((strain_principal(2,2)+abs(strain_principal(2,2)))/2)^2);   
end