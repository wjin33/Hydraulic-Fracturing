% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function  stiffnessMatrix(PROP,LTAN,NLTAN,UPDATE)
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

global CONNEC NODES PSI XYZ STATEV DISPDD DISPTD FORCE GKF connec_frac xyz_frac statev_frac LocalElement

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

            nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
            [NLEquivStrain,scale] = computeNonlocalEquivalentStrain( nonlocal_table );
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;

            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal);
            
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
            
            % Update residual Force by minus internal forces
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];

            %Update tangent/secant stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                        % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;        % Sum up all Gauss point contribution/% Assemble the global stiffness
                
                if NLTAN
                    [dC_domega1, dC_domega2] = Gauss_tan_derivative(PROP,damageN);
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    
                    if ( lcoeff(1) == 0 && lcoeff(2) == 0 )
                        continue;
                    end
                    
                    NNL = size(nonlocal_table,1);
                    for inl = 1:NNL
                        
                        iele_nl = nonlocal_table(inl,1);
                        igp_nl = nonlocal_table(inl,2);
                        weight_nl = nonlocal_table(inl,3);
                        volume_nl = STATEV{iele_nl}{igp_nl}.volume;
                        
                        N1_nl  = CONNEC(iele_nl,2);                                                  % Node 1 for the element that has gauss points is located inside the nonlcal influence zone
                        N2_nl  = CONNEC(iele_nl,3);                                                  % Node 2 for the element that has gauss points is located inside the nonlcal influence zone
                        N3_nl  = CONNEC(iele_nl,4);                                                  % Node 3 for the element that has gauss points is located inside the nonlcal influence zone
                        N4_nl  = CONNEC(iele_nl,5);                                                  % Node 4 for the element that has gauss points is located inside the nonlcal influence zone
                        NN_nl  = NODES([N1_nl N2_nl N3_nl N4_nl]',:);                                          % Nodal data for current element             

                        local_nl  = [N1_nl*2-1 N1_nl*2 N2_nl*2-1 N2_nl*2 N3_nl*2-1 N3_nl*2 N4_nl*2-1 N4_nl*2];             % Traditional index locations
                        xyz_nl = [XYZ(N1_nl,2) XYZ(N1_nl,3);
                                  XYZ(N2_nl,2) XYZ(N2_nl,3);
                                  XYZ(N3_nl,2) XYZ(N3_nl,3);
                                  XYZ(N4_nl,2) XYZ(N4_nl,3);];                                % Nodal coordinate matrix for the elements that the nonlocal gauss point is located inside

                        gp_nl = STATEV{iele_nl}{igp_nl}.natural_coodinates;                 % Gauss points
%                         W_nl = STATEV{iele_nl}{igp_nl}.gauss_weight;                              % Gauss weights               
                        xi_nl = gp_nl(1,1); eta_nl = gp_nl(2,1);                                % Gauss points
                        [Nxy_nl, ~] = Shape_Function(xi_nl, eta_nl, xyz_nl);                % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian
                        
                        N_nl  = 1/4*[(1-xi_nl)*(1-eta_nl);(1+xi_nl)*(1-eta_nl);...                     % Shape functions
                                     (1+xi_nl)*(1+eta_nl);(1-xi_nl)*(1+eta_nl)];
                        
                        Bu_nl = [Nxy_nl(1,1)   0             Nxy_nl(1,2)   0             Nxy_nl(1,3)   0             Nxy_nl(1,4)   0;...
                                 0             Nxy_nl(2,1)   0             Nxy_nl(2,2)   0             Nxy_nl(2,3)   0             Nxy_nl(2,4);...
                                 Nxy_nl(2,1)   Nxy_nl(1,1)   Nxy_nl(2,2)   Nxy_nl(1,2)   Nxy_nl(2,3)   Nxy_nl(1,3)   Nxy_nl(2,4)   Nxy_nl(1,4)];
                        Benr_nl = [];
                        
                        index = 1;
                        iLoc_nl = 9 ;
                        for iN = 1:4
                            if NN_nl(iN,2) ~= 0
                                psi1 = PSI(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI(N4_nl);                                              % Psi level set value at node 4
                                psi  = N_nl(1)*psi1+N_nl(2)*psi2+N_nl(3)*psi3+N_nl(4)*psi4;     % Psi level set value at current gauss point
    
                                Hgp = sign(psi);                                        % Heaviside value at current nonlocal gauss point
                                Hi  = NN_nl(iN,3);                                         % Nodal Heaviside value
                                H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

                                Ba = [Nxy_nl(1,iN)*H    0;
                                      0                 Nxy_nl(2,iN)*H;
                                      Nxy_nl(2,iN)*H    Nxy_nl(1,iN)*H];
                                  
                                Benr_nl(:,(2*index-1):(2*index)) = Ba;
                                index = index+1;
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [2*NN_nl(iN,2)-1 2*NN_nl(iN,2)];
                                iLoc_nl = iLoc_nl+2;
                            end
                        end
                        
              
                        DSPD_nl=DISPDD(local_nl);
                        B_nl = [Bu_nl Benr_nl];
                        deps_nl = B_nl*DSPD_nl;
                        DEPS_nl =[deps_nl(1,1); deps_nl(2,1); 0; 0.5*deps_nl(3,1);];
                        strain_nl = STATEV{iele_nl}{igp_nl}.strain;
                        strain_nl = strain_nl + DEPS_nl ;
                        
                        equivalentEPS_nl = STATEV{iele_nl}{igp_nl}.EquivStrain; 
                          
                        eqeps_1t = PROP.eqeps_1t;
                        eqeps_2t = PROP.eqeps_2t;
                        eqeps_1s = PROP.eqeps_1s;
                        
                        if (equivalentEPS_nl(1) ~= 0 && lcoeff(1) ~= 0)
                            dequieps_deps = [strain_nl(1)  0             strain_nl(4)*(eqeps_1t/eqeps_1s)^2]./equivalentEPS_nl(1);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(1)*detJ*W*Bu'*dC_domega1*CBu_nl;
                        end
                        if (equivalentEPS_nl(2) ~= 0 && lcoeff(2) ~= 0 )
                            dequieps_deps = [0             strain_nl(2)  strain_nl(4)*(eqeps_2t/eqeps_1s)^2]./equivalentEPS_nl(2);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(2)*detJ*W*Bu'*dC_domega2*CBu_nl;
                        end
 
                    end
                    
                end
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
            
            nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
            [NLEquivStrain,scale] = computeNonlocalEquivalentStrain( nonlocal_table );
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal);
                        
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
            
            % Update residual FORCE by minus internal FORCEs
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];

            % Tangent stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                     % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;                 % Sum up all Gauss point contribution                
                
                if NLTAN
                    [dC_domega1, dC_domega2] = Gauss_tan_derivative(PROP,damageN);
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    
                    if ( lcoeff(1) == 0 && lcoeff(2) == 0 )
                        continue;
                    end
                    
                    NNL = size(nonlocal_table,1);
                    for inl = 1:NNL
                        
                        iele_nl = nonlocal_table(inl,1);
                        igp_nl = nonlocal_table(inl,2);
                        weight_nl = nonlocal_table(inl,3);
                        volume_nl = STATEV{iele_nl}{igp_nl}.volume;
                        
                        N1_nl  = CONNEC(iele_nl,2);                                                  % Node 1 for the element that has gauss points is located inside the nonlcal influence zone
                        N2_nl  = CONNEC(iele_nl,3);                                                  % Node 2 for the element that has gauss points is located inside the nonlcal influence zone
                        N3_nl  = CONNEC(iele_nl,4);                                                  % Node 3 for the element that has gauss points is located inside the nonlcal influence zone
                        N4_nl  = CONNEC(iele_nl,5);                                                  % Node 4 for the element that has gauss points is located inside the nonlcal influence zone
                        NN_nl  = NODES([N1_nl N2_nl N3_nl N4_nl]',:);                                          % Nodal data for current element

                        local_nl  = [N1_nl*2-1 N1_nl*2   N2_nl*2-1 N2_nl*2   N3_nl*2-1 N3_nl*2   N4_nl*2-1 N4_nl*2];             % Traditional index locations
                                    
                        xyz_nl = [XYZ(N1_nl,2) XYZ(N1_nl,3);
                                  XYZ(N2_nl,2) XYZ(N2_nl,3);
                                  XYZ(N3_nl,2) XYZ(N3_nl,3);
                                  XYZ(N4_nl,2) XYZ(N4_nl,3);];                                % Nodal coordinate matrix for the elements that the nonlocal gauss point is located inside

                        gp_nl = STATEV{iele_nl}{igp_nl}.natural_coodinates;                 % Gauss points
                        xi_nl = gp_nl(1,1); eta_nl = gp_nl(2,1);                                % Gauss points
                        [Nxy_nl, ~] = Shape_Function(xi_nl, eta_nl, xyz_nl);                % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian  
                        
                        N_nl  = 1/4*[(1-xi_nl)*(1-eta_nl);(1+xi_nl)*(1-eta_nl);...                     % Shape functions
                                  (1+xi_nl)*(1+eta_nl);(1-xi_nl)*(1+eta_nl)];

                        Bu_nl = [Nxy_nl(1,1)   0             Nxy_nl(1,2)   0             Nxy_nl(1,3)   0             Nxy_nl(1,4)   0;...
                                 0             Nxy_nl(2,1)   0             Nxy_nl(2,2)   0             Nxy_nl(2,3)   0             Nxy_nl(2,4);...
                                 Nxy_nl(2,1)   Nxy_nl(1,1)   Nxy_nl(2,2)   Nxy_nl(1,2)   Nxy_nl(2,3)   Nxy_nl(1,3)   Nxy_nl(2,4)   Nxy_nl(1,4)];
                        Benr_nl = [];     
                        index = 1;
                        iLoc_nl = 9 ;
                        for iN = 1:4
                            if NN_nl(iN,2) ~= 0
                                psi1 = PSI(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI(N4_nl);                                              % Psi level set value at node 4
                                psi  = N_nl(1)*psi1+N_nl(2)*psi2+N_nl(3)*psi3+N_nl(4)*psi4;     % Psi level set value at current gauss point
    
                                Hgp = sign(psi);                                        % Heaviside value at current nonlocal gauss point
                                Hi  = NN_nl(iN,3);                                         % Nodal Heaviside value
                                H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

%                                 Hgp = abs(psi);                                        % Heaviside value at current gauss point
%                                 Hi  = abs(PSI(NN(iN,1)));                                         % Nodal Heaviside value
%                                 H   = (Hgp-Hi); 
                    
                                Ba = [Nxy_nl(1,iN)*H    0;
                                      0                 Nxy_nl(2,iN)*H;
                                      Nxy_nl(2,iN)*H    Nxy_nl(1,iN)*H];
                                  
                                Benr_nl(:,(2*index-1):(2*index)) = Ba;
                                index = index+1;
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [2*NN_nl(iN,2)-1 2*NN_nl(iN,2)];
                                iLoc_nl = iLoc_nl+2;
                            end
                        end
                        
                        B_nl = [Bu_nl Benr_nl]; 
                        DSPD_nl=DISPDD(local_nl);
                        deps_nl = B_nl*DSPD_nl;
                        DEPS_nl =[deps_nl(1,1); deps_nl(2,1); 0; 0.5*deps_nl(3,1);];
                        strain_nl = STATEV{iele_nl}{igp_nl}.strain;
                        strain_nl = strain_nl + DEPS_nl ;
                        
                        equivalentEPS_nl = STATEV{iele_nl}{igp_nl}.EquivStrain; 
                          
                        eqeps_1t = PROP.eqeps_1t;
                        eqeps_2t = PROP.eqeps_1t;
                        eqeps_1s = PROP.eqeps_1t;
                        
                        if (equivalentEPS_nl(1) ~= 0 && lcoeff(1) ~= 0)
                            dequieps_deps = [strain_nl(1)  0             strain_nl(4)*(eqeps_1t/eqeps_1s)^2]./equivalentEPS_nl(1);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(1)*detJ*W*Bu'*dC_domega1*CBu_nl;
                        end
                        if (equivalentEPS_nl(2) ~= 0 && lcoeff(2) ~= 0 )
                            dequieps_deps = [0             strain_nl(2)  strain_nl(4)*(eqeps_2t/eqeps_1s)^2]./equivalentEPS_nl(2);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(2)*detJ*W*Bu'*dC_domega2*CBu_nl;
                        end
                        
                    end
                    
                end
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
            
            nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
            [NLEquivStrain,scale] = computeNonlocalEquivalentStrain( nonlocal_table );
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,isLocal);
                        
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
            
            % Update residual Force by minus internal forces
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];
            
            % Tangent stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                     % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;                 % Sum up all Gauss point contribution                
                
                if NLTAN
                    [dC_domega1, dC_domega2] = Gauss_tan_derivative(PROP,damageN);
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    
                    if ( lcoeff(1) == 0 && lcoeff(2) == 0 )
                        continue;
                    end
                    
                    NNL = size(nonlocal_table,1);
                    for inl = 1:NNL
                        
                        iele_nl = nonlocal_table(inl,1);
                        igp_nl = nonlocal_table(inl,2);
                        weight_nl = nonlocal_table(inl,3);
                        volume_nl = STATEV{iele_nl}{igp_nl}.volume;
                        
                        N1_nl  = CONNEC(iele_nl,2);                                                  % Node 1 for the element that has gauss points is located inside the nonlcal influence zone
                        N2_nl  = CONNEC(iele_nl,3);                                                  % Node 2 for the element that has gauss points is located inside the nonlcal influence zone
                        N3_nl  = CONNEC(iele_nl,4);                                                  % Node 3 for the element that has gauss points is located inside the nonlcal influence zone
                        N4_nl  = CONNEC(iele_nl,5);                                                  % Node 4 for the element that has gauss points is located inside the nonlcal influence zone
                        NN_nl  = NODES([N1_nl N2_nl N3_nl N4_nl]',:);                                          % Nodal data for current element

                        local_nl  = [N1_nl*2-1 N1_nl*2   N2_nl*2-1 N2_nl*2   N3_nl*2-1 N3_nl*2   N4_nl*2-1 N4_nl*2];             % Traditional index locations
                        xyz_nl = [XYZ(N1_nl,2) XYZ(N1_nl,3);
                                  XYZ(N2_nl,2) XYZ(N2_nl,3);
                                  XYZ(N3_nl,2) XYZ(N3_nl,3);
                                  XYZ(N4_nl,2) XYZ(N4_nl,3);];                                % Nodal coordinate matrix for the elements that the nonlocal gauss point is located inside

                        gp_nl = STATEV{iele_nl}{igp_nl}.natural_coodinates;                 % Gauss points
                        xi_nl = gp_nl(1,1); eta_nl = gp_nl(2,1);                                % Gauss points
                        [Nxy_nl, ~] = Shape_Function(xi_nl, eta_nl, xyz_nl);                % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian  
                        
                        N_nl  = 1/4*[(1-xi_nl)*(1-eta_nl);(1+xi_nl)*(1-eta_nl);...                     % Shape functions
                                  (1+xi_nl)*(1+eta_nl);(1-xi_nl)*(1+eta_nl)];

                        Bu_nl = [Nxy_nl(1,1)   0             Nxy_nl(1,2)   0             Nxy_nl(1,3)   0             Nxy_nl(1,4)   0;...
                                 0             Nxy_nl(2,1)   0             Nxy_nl(2,2)   0             Nxy_nl(2,3)   0             Nxy_nl(2,4);...
                                 Nxy_nl(2,1)   Nxy_nl(1,1)   Nxy_nl(2,2)   Nxy_nl(1,2)   Nxy_nl(2,3)   Nxy_nl(1,3)   Nxy_nl(2,4)   Nxy_nl(1,4)];
                        Benr_nl = [];     
                        index = 1;
                        iLoc_nl = 9 ;
                        for iN = 1:4
                            if NN_nl(iN,2) ~= 0
                                psi1 = PSI(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI(N4_nl);                                              % Psi level set value at node 4
                                psi  = N_nl(1)*psi1+N_nl(2)*psi2+N_nl(3)*psi3+N_nl(4)*psi4;     % Psi level set value at current gauss point
    
                                Hgp = sign(psi);                                        % Heaviside value at current nonlocal gauss point
                                Hi  = NN_nl(iN,3);                                         % Nodal Heaviside value
                                H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

                                Ba = [Nxy_nl(1,iN)*H    0;
                                      0                 Nxy_nl(2,iN)*H;
                                      Nxy_nl(2,iN)*H    Nxy_nl(1,iN)*H];
                                  
                                Benr_nl(:,(2*index-1):(2*index)) = Ba;
                                index = index+1;
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [2*NN_nl(iN,2)-1 2*NN_nl(iN,2)-0];
                                iLoc_nl = iLoc_nl+2;
                            end
                        end
                        
                        B_nl = [Bu_nl Benr_nl]; 
                        DSPD_nl=DISPDD(local_nl);
                        deps_nl = B_nl*DSPD_nl;
                        DEPS_nl =[deps_nl(1,1); deps_nl(2,1); 0; 0.5*deps_nl(3,1);];
                        strain_nl = STATEV{iele_nl}{igp_nl}.strain;
                        strain_nl = strain_nl + DEPS_nl ;
                        
                        equivalentEPS_nl = STATEV{iele_nl}{igp_nl}.EquivStrain; 
                          
                        eqeps_1t = PROP.eqeps_1t;
                        eqeps_2t = PROP.eqeps_1t;
                        eqeps_1s = PROP.eqeps_1t;
                        
                        if (equivalentEPS_nl(1) ~= 0 && lcoeff(1) ~= 0)
                            dequieps_deps = [strain_nl(1)  0             strain_nl(4)*(eqeps_1t/eqeps_1s)^2]./equivalentEPS_nl(1);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(1)*detJ*W*Bu'*dC_domega1*CBu_nl;
                        end
                        if (equivalentEPS_nl(2) ~= 0 && lcoeff(2) ~= 0 )
                            dequieps_deps = [0             strain_nl(2)  strain_nl(4)*(eqeps_2t/eqeps_1s)^2]./equivalentEPS_nl(2);
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff(2)*detJ*W*Bu'*dC_domega2*CBu_nl;
                        end
                    end
                end
            end
        end
        
        I_frac = find( connec_frac(:,1) == iElem );
        
        N1_frac  = connec_frac(I_frac,2);                                              % Node 1 for current fracture segment
        N2_frac  = connec_frac(I_frac,3);                                              % Node 2 for current fracture segment
        
        X1_frac = xyz_frac(N1_frac,4);
        X2_frac = xyz_frac(N2_frac,4);
        
        tan_frac = xyz_frac(N2_frac,2:3) - xyz_frac(N1_frac,2:3);
        tan_frac = tan_frac./sqrt(sum(tan_frac.^2));
        norm_frac = [-tan_frac(2) tan_frac(1)];
        
        lambda = [tan_frac; norm_frac];
        
        localD_enr  = localD(1,9:16);
        
        Ngp = size(statev_frac{I_frac},2);
        
        for ig = 1:Ngp

            gp_ele = statev_frac{I_frac}{ig}.natural_coordinate_ele;            % Gauss points
            W = statev_frac{I_frac}{ig}.gauss_weight;                        % Gauss weights
            histo_jump = statev_frac{I_frac}{ig}.jump; 
            
            xi = gp_ele(1,1); eta = gp_ele(2,1);                                % Gauss points
            
            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                     % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];
                  
            detJ_frac = det([-0.5 0.5]*[X1_frac; X2_frac]);
 

            NT = [N(1)   0     N(2)  0     N(3)  0     N(4)  0 ;
                  0      N(1)  0     N(2)  0     N(3)  0     N(4);];
     
            jump = lambda*NT*DISPTD(localD_enr);
            fracwidth = jump(2);
            
            Traction =  cohesiveLaw(jump, histo_jump);
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

                FORCE(localD_enr) = FORCE(localD_enr) - W*detJ_frac*NT'*lambda'*Traction;
                    
                if LTAN
                    M_czm =  cohesiveStiffness(jump, histo_jump);
                    K_CZM = W*NT'*lambda'*M_czm*lambda*NT*detJ_frac;  %% cohesive zone section
                    
                    GKF(localD_enr,localD_enr) = GKF(localD_enr,localD_enr) + K_CZM;
                end
        end
    end 
end

end

function [stressN, DamageN, kappaN]=Gauss_sig_update(PROP,EquivStrain,strainN,kappa,Damage,isLocal)
% Inputs:
% PROP = [E0, nu0, epsilon_0, epsilon_f];
% D = 4*4 elastic stiffness matrix
% stress = [s11, s22, s12,   s33];
% strain = [e11, e22, 2*e12, e33];
% Plane strain problem e33=0
%%
    E11  = PROP.E11;
    E22  = PROP.E22;
    nu12 = PROP.nu12;
    nu21 = nu12*E11/E22;
    nu23 = PROP.nu23;
    G12 = PROP.G12;

    eqeps_1t = PROP.eqeps_1t;
    eqeps_2t = PROP.eqeps_2t;
 
    alpha_1t = PROP.alpha_1t;
    alpha_2t = PROP.alpha_2t;
    if ~isLocal
        if ( EquivStrain(1) > kappa(1) )
            Damage(1) = 1-exp(-(EquivStrain(1) - eqeps_1t)/alpha_1t);
            kappa(1) = EquivStrain(1);
        end
        if ( EquivStrain(2) > kappa(2) )
            Damage(2) = 1-exp(-(EquivStrain(2) - eqeps_2t)/alpha_2t);
            kappa(2) = EquivStrain(2);
        end
    end
    
    DamageN = Damage;
    kappaN = kappa;
    
    omega1 = DamageN(1);
    omega2 = DamageN(2);

%     omega1 = 0;
%     omega2 = 0;
    
    
    MATC=zeros(4,4);

    nu21=E22*nu12/E11;
    D= (1-omega2)*nu23^2+2*(1-omega1)*(1-omega2)*nu12*nu21*nu23+(1-omega1)*(2-omega2)*nu12*nu21-1;
    MATC(1,1)=E11*(1-omega1)*((1-omega2)*nu23^2-1)/D;
    MATC(1,2)=-E11*nu21*(1-omega1)*(1-omega2)*(1+nu23)/D;
    MATC(1,3)=-E11*nu21*(1-omega1)*(1+(1-omega2)*nu23)/D;
    MATC(1,4)=0;
    MATC(2,1)=MATC(1,2);
    MATC(2,2)=E22*(1-omega2)*((1-omega1)*nu12*nu21-1)/D;
    MATC(2,3)=-E22*(1-omega2)*(nu23+(1-omega1)*nu12*nu21)/D;
    MATC(2,4)=0;
    MATC(3,1)=MATC(1,3);
    MATC(3,2)=MATC(2,3);
    MATC(3,3)=E22*(1-omega2)*(1-omega1)*(nu12*nu21-1)/D;
    MATC(3,4)=0;
    MATC(4,1)=0;
    MATC(4,2)=0;
    MATC(4,3)=0;
    MATC(4,4)=G12*(1-omega1)*(1-omega2);
    
%     alpha = BiotCoefficient_tan_update(PROP);

    stressN = MATC*[strainN(1:3,1); 2*strainN(4,1)];    %updated stress

end

function [nonlocal_equ_eps,scale] = computeNonlocalEquivalentStrain( nonlocal_table )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global  STATEV

    nonlocal_equ_eps = [0; 0;];
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

function [MATC]=Gauss_tan_update(PROP,Damage)
% Inputs:
% Plane strain problem e33=0
%%
    E11  = PROP.E11;
    E22  = PROP.E22;
    nu12 = PROP.nu12;
    nu21 = nu12*E11/E22;
    nu23 = PROP.nu23;
    G12 = PROP.G12;
    
    omega1 = Damage(1);
    omega2 = Damage(2);

%     omega1 = 0;
%     omega2 = 0;

    MATC=zeros(3,3);

    nu21=E22*nu12/E11;
    D= (1-omega2)*nu23^2+2*(1-omega1)*(1-omega2)*nu12*nu21*nu23+(1-omega1)*(2-omega2)*nu12*nu21-1;
    MATC(1,1)=E11*(1-omega1)*((1-omega2)*nu23^2-1)/D;
    MATC(1,2)=-E11*nu21*(1-omega1)*(1-omega2)*(1+nu23)/D;
    MATC(1,3)=0;
    MATC(2,1)=MATC(1,2);
    MATC(2,2)=E22*(1-omega2)*((1-omega1)*nu12*nu21-1)/D;
    MATC(2,3)=0;
    MATC(3,1)=0;
    MATC(3,2)=0;
    MATC(3,3)=G12*(1-omega1)*(1-omega2);
end

function [dC_domega1, dC_domega2]=Gauss_tan_derivative(PROP,Damage)
% Inputs:
% Plane strain problem e33=0
%%

    E11  = PROP.E11;
    E22  = PROP.E22;
    nu12 = PROP.nu12;
    nu21 = nu12*E11/E22;
    nu23 = PROP.nu23;
    G12 = PROP.G12;
    
    omega1 = Damage(1);
    omega2 = Damage(2);
    
    nu21=E22*nu12/E11;
    
    D= (1-omega2)*nu23^2+2*(1-omega1)*(1-omega2)*nu12*nu21*nu23+(1-omega1)*(2-omega2)*nu12*nu21-1;
    C11=E11*(1-omega1)*((1-omega2)*nu23^2-1);
    C12=-E11*nu21*(1-omega1)*(1-omega2)*(1+nu23);
    C22=E22*(1-omega2)*((1-omega1)*nu12*nu21-1);
%     C33=G12*(1-omega1)*(1-omega2);
        
    dD_domega1 = -2*(1-omega2)*nu12*nu21*nu23 - (2-omega2)*nu12*nu21;
    dD_domega2 = -nu23^2 - 2*(1-omega1)*nu12*nu21*nu23 - (1-omega1)*nu12*nu21;
    
    dC11_domega1 = - E11*((1-omega2)*nu23^2-1);
    dC11_domega2 = - E11*(1-omega1)*nu23^2;
    
    dC22_domega1 = - E22*(1-omega2)*nu12*nu21;
    dC22_domega2 = - E22*((1-omega1)*nu12*nu21-1);
    
    dC33_domega1 = - G12*(1-omega2);
    dC33_domega2 = - G12*(1-omega1);
    
    dC12_domega1 = E11*nu21*(1-omega2)*(1+nu23);
    dC12_domega2 = E11*nu21*(1-omega1)*(1+nu23);
    
    dC_domega1 = zeros(3,3);
    dC_domega1(1,1)=( dC11_domega1*D - dD_domega1*C11 )/D^2;
    dC_domega1(1,2)=( dC12_domega1*D - dD_domega1*C12 )/D^2;
    dC_domega1(1,3)=0;
    dC_domega1(2,1)=dC_domega1(1,2);
    dC_domega1(2,2)=( dC22_domega1*D - dD_domega1*C22 )/D^2;
    dC_domega1(2,3)=0;
    dC_domega1(3,1)=0;
    dC_domega1(3,2)=0;
    dC_domega1(3,3)=dC33_domega1;
    
    dC_domega2 = zeros(3,3);
    dC_domega2(1,1)=( dC11_domega2*D - dD_domega2*C11 )/D^2;
    dC_domega2(1,2)=( dC12_domega2*D - dD_domega2*C12 )/D^2;
    dC_domega2(1,3)=0;
    dC_domega2(2,1)=dC_domega2(1,2);
    dC_domega2(2,2)=( dC22_domega2*D - dD_domega2*C22 )/D^2;
    dC_domega2(2,3)=0;
    dC_domega2(3,1)=0;
    dC_domega2(3,2)=0;
    dC_domega2(3,3)=dC33_domega2;
end

function [lcoeff] = LocalCoefficient(PROP,kappa,NLEquivStrain)

    eqeps_1t = PROP.eqeps_1t;
    eqeps_2t = PROP.eqeps_2t;
 
    alpha_1t = PROP.alpha_1t;
    alpha_2t = PROP.alpha_2t;
    
    lcoeff =  zeros(2,1);
    
    if (NLEquivStrain(1) == kappa(1))
        lcoeff(1) = exp(-(NLEquivStrain(1) - eqeps_1t)/alpha_1t)/alpha_1t;
    end
    if (NLEquivStrain(2) == kappa(2))
        lcoeff(2) = exp(-(NLEquivStrain(2) - eqeps_2t)/alpha_2t)/alpha_2t;
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








