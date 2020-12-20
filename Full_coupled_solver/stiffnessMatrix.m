% Written By: Wencheng Jin, Idaho National Laboratory (2020)
% Website: https://sites.google.com/view/wenchengjin/software
% Email: wencheng.jin@inl.gov

function  stiffnessMatrix(PROP,LTAN,NLTAN,UPDATE,Theta,deltaT)
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

global CONNEC NODES PSI XYZ STATEV DISPDD DISPTD PREFTD FORCE GKF  connec_frac xyz_frac statev_frac

nCrack = size(connec_frac,1);

for iElem = 1:size(CONNEC,1)
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    
    localD  = [N1*3-2 N1*3-1 N2*3-2 N2*3-1 N3*3-2 N3*3-1 N4*3-2 N4*3-1];    % Traditional index locations for displacement
    localP  = [N1*3 N2*3 N3*3 N4*3];                                        % Traditional index locations for pressure
    iLocD   = 9;                                                            % Next index location
    iLocP   = 5;
    
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
%             N  = 1/4*[(1+xi)*(1+eta) (1-xi)*(1+eta) ...                    % Shape functions
%                       (1-xi)*(1-eta) (1+xi)*(1-eta)];

            [Nxy, detJ] = Shape_Function(xi, eta, xyz);                     % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian             
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];
              
            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;
            
            ELEPORP=DISPTD(localP);
            K = HydroConductivity_tan_update(PROP);
            velocityN = -K*Nxy*ELEPORP;
            GPporepressure = N*ELEPORP;
            
            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,scale] = computeNonlocalEquivalentStrain( nonlocal_table );
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;

            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,GPporepressure);
            
            % Update plastic variables once convergence is reached within
            % one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
                STATEV{iElem}{i}.fluidVelocity = velocityN;
                continue;
            end
           
            % Update residual Force by minus internal forces
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];
            
            M_coef=BiotModulus_tan_update(PROP);
            H_coef=HydroConductivity_tan_update(PROP);
            Q_coef=BiotCoefficient_tan_update(PROP);
                
            M = W*N'*M_coef*N*detJ;
            H = W*Nxy'*H_coef*Nxy*detJ;
            Q = W*Bu'*Q_coef'*N*detJ;
            QT = W*N'*Q_coef*Bu*detJ;
%             
%             FORCE(localP) = FORCE(localP)+QT*( 1/(Theta*deltaT)*PREFTD(localD)+(1/Theta-1)*DISPTD_Rate(localD) )...
%                                           +M*( 1/(Theta*deltaT)*PREFTD(localP)+(1/Theta-1)*DISPTD_Rate(localP) )...
%                                           -1/(Theta*deltaT)*QT*DISPTD(localD) - (1/(Theta*deltaT)*M+H)*DISPTD(localP);
            FORCE(localP) = FORCE(localP) - QT*(DISPTD(localD)-PREFTD(localD)) - M*(DISPTD(localP)-PREFTD(localP))...
                            -deltaT*H*(Theta*DISPTD(localP)+(1-Theta)*PREFTD(localP)) ;

            %Update tangent/secant stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                        % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;        % Sum up all Gauss point contribution/% Assemble the global stiffness
%                 
%                 GKF(localP,localP) = GKF(localP,localP) + M/(Theta*deltaT) + H;
%                 GKF(localD,localP) = GKF(localD,localP) - Q;
%                 GKF(localP,localD) = GKF(localP,localD) + 1/(Theta*deltaT)*QT;
                GKF(localP,localP) = GKF(localP,localP) + M + H*(Theta*deltaT);
                GKF(localD,localP) = GKF(localD,localP) - Q;
                GKF(localP,localD) = GKF(localP,localD) + QT;

                if NLTAN
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    if ( lcoeff == 0 )
                        continue;
                    end
                    [dC_domega] = Gauss_tan_derivative(PROP);
                    
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

                        local_nl  = [N1_nl*3-2 N1_nl*3-1 N2_nl*3-2 N2_nl*3-1 N3_nl*3-2 N3_nl*3-1 N4_nl*3-2 N4_nl*3-1];             % Traditional index locations
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
                                icrack =  NN_nl(iN,4);
                                psi1 = PSI{icrack}(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI{icrack}(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI{icrack}(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI{icrack}(N4_nl);                                              % Psi level set value at node 4
                                psi  = N_nl(1)*psi1+N_nl(2)*psi2+N_nl(3)*psi3+N_nl(4)*psi4;     % Psi level set value at current gauss point
    
                                Hgp = sign(psi);                                        % Heaviside value at current nonlocal gauss point
                                Hi  = NN_nl(iN,3);                                         % Nodal Heaviside value
                                H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

                                Ba = [Nxy_nl(1,iN)*H    0;
                                      0                 Nxy_nl(2,iN)*H;
                                      Nxy_nl(2,iN)*H    Nxy_nl(1,iN)*H];
                                  
                                Benr_nl(:,(2*index-1):(2*index)) = Ba;
                                index = index+1;
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [3*NN_nl(iN,2)-2 3*NN_nl(iN,2)-1];
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
                        
                        if (equivalentEPS_nl ~= 0 && lcoeff ~= 0)
                            [esp_pri_dir,esp_pri_val]= eigs([strain_nl(1,1) strain_nl(4,1); strain_nl(4,1) strain_nl(2,1)]);
                            angle_eps_pri=[0;0];
                            if (esp_pri_val(1,1) > 0)
                                angle_eps_pri(1) = esp_pri_val(1,1);
                            end
                            if (esp_pri_val(2,2) > 0)
                                angle_eps_pri(2) = esp_pri_val(2,2);
                            end
                            dequieps_deps = [esp_pri_dir(1,1)*angle_eps_pri(1)   esp_pri_dir(2,2)*angle_eps_pri(2)  esp_pri_dir(1,2)*angle_eps_pri(2)+esp_pri_dir(2,1)*angle_eps_pri(1)]./equivalentEPS_nl;
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff*detJ*W*Bu'*dC_domega*CBu_nl;
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
            Nenr = [];
            
            Benru = [];
            Benrp = [];
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            index = 1;
            R_weight = 0;
            for iN = 1:4
                if NN(iN,2) ~= 0
                    icrack =  NN(iN,4);
                    psi1 = PSI{icrack}(N1);                                         % Psi level set value at node 1
                    psi2 = PSI{icrack}(N2);                                         % Psi level set value at node 2
                    psi3 = PSI{icrack}(N3);                                         % Psi level set value at node 3
                    psi4 = PSI{icrack}(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
    
                    Hgp = sign(psi);                                        % Heaviside value at current gauss point
                    Hi  = NN(iN,3);                                         % Nodal Heaviside value
                    H   = (Hgp-Hi)/2;                                       % Shifted Heaviside value
                    
                    R_weight = R_weight + N(iN);
                    
                    Hgp = abs(psi);                                        % Heaviside value at current gauss point
                    Hi  = abs(PSI{icrack}(NN(iN,1)));                              % Nodal Heaviside value
                    D   = (Hgp-Hi);                                        % Shifted distance value

                    Benru(:,(2*index-1):(2*index)) = [Nxy(1,iN)*H    0;
                                                      0              Nxy(2,iN)*H;
                                                      Nxy(2,iN)*H    Nxy(1,iN)*H];
                    Benrp(:,index) = [Nxy(1,iN)*D;Nxy(2,iN)*D];
                    Nenr(:,index) = N(1,iN)*D;
                    index = index+1;
                    if (i == 1)
                        localD(iLocD:(iLocD+1)) = [3*NN(iN,2)-2 3*NN(iN,2)-1];
                        iLocD = iLocD+2;
                        localP(iLocP) = 3*NN(iN,2);
                        iLocP = iLocP+1;
                    end
                end    
            end

            Bu = [Bu Benru];
            Bp = [Nxy Benrp.*R_weight];
            Nen = [N Nenr.*R_weight];
            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;
            
            ELEPORP=DISPTD(localP);
            K = HydroConductivity_tan_update(PROP);
            velocityN = -K*Bp*ELEPORP;
            GPporepressure = Nen*ELEPORP;
            
            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,scale] = computeNonlocalEquivalentStrain( nonlocal_table );
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,GPporepressure);
                        
            % Update plastic variables once convergence is reached within one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
                STATEV{iElem}{i}.fluidVelocity = velocityN;
               continue;
            end
            
            % Update residual Force by minus internal forces
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];
            
            M_coef=BiotModulus_tan_update(PROP);
            H_coef=HydroConductivity_tan_update(PROP);
            Q_coef=BiotCoefficient_tan_update(PROP);
                
            M = W*Nen'*M_coef*Nen*detJ;
            H = W*Bp'*H_coef*Bp*detJ;
            Q = W*Bu'*Q_coef'*Nen*detJ;
            QT = W*Nen'*Q_coef*Bu*detJ;
            
            FORCE(localP) = FORCE(localP) - QT*(DISPTD(localD)-PREFTD(localD)) - M*(DISPTD(localP)-PREFTD(localP))...
                            -deltaT*H*(Theta*DISPTD(localP)+(1-Theta)*PREFTD(localP)) ;
            
%             FORCE(localP) = FORCE(localP)+QT*( 1/(Theta*deltaT)*PREFTD(localD)+(1/Theta-1)*DISPTD_Rate(localD) )...
%                                           +M*( 1/(Theta*deltaT)*PREFTD(localP)+(1/Theta-1)*DISPTD_Rate(localP) )...
%                                           -1/(Theta*deltaT)*QT*DISPTD(localD) - (1/(Theta*deltaT)*M+H)*DISPTD(localP);

            % Tangent stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                     % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;                 % Sum up all Gauss point contribution                
%                 GKF(localP,localP) = GKF(localP,localP) + 1/(Theta*deltaT)*M + H;                
%                 GKF(localD,localP) = GKF(localD,localP) - Q;
%                 GKF(localP,localD) = GKF(localP,localD) + 1/(Theta*deltaT)*QT;
                GKF(localP,localP) = GKF(localP,localP) + M + H*(Theta*deltaT);
                GKF(localD,localP) = GKF(localD,localP) - Q;
                GKF(localP,localD) = GKF(localP,localD) + QT;
                
                if NLTAN
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    if ( lcoeff == 0 )
                        continue;
                    end
                    [dC_domega] = Gauss_tan_derivative(PROP);
                    
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

                        local_nl  = [N1_nl*3-2 N1_nl*3-1 N2_nl*3-2 N2_nl*3-1 N3_nl*3-2 N3_nl*3-1 N4_nl*3-2 N4_nl*3-1];             % Traditional index locations
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
                                icrack =  NN_nl(iN,4);
                                psi1 = PSI{icrack}(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI{icrack}(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI{icrack}(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI{icrack}(N4_nl); 
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
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [3*NN_nl(iN,2)-2 3*NN_nl(iN,2)-1];
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
                          
                        if (equivalentEPS_nl ~= 0 && lcoeff ~= 0)
                            [esp_pri_dir,esp_pri_val]= eigs([strain_nl(1,1) strain_nl(4,1); strain_nl(4,1) strain_nl(2,1)]);
                            angle_eps_pri=[0;0];
                            if (esp_pri_val(1,1) > 0)
                                angle_eps_pri(1) = esp_pri_val(1,1);
                            end
                            if (esp_pri_val(2,2) > 0)
                                angle_eps_pri(2) = esp_pri_val(2,2);
                            end
                            dequieps_deps = [esp_pri_dir(1,1)*angle_eps_pri(1)   esp_pri_dir(2,2)*angle_eps_pri(2)  esp_pri_dir(1,2)*angle_eps_pri(2)+esp_pri_dir(2,1)*angle_eps_pri(1)]./equivalentEPS_nl;
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff*detJ*W*Bu'*dC_domega*CBu_nl;
                        end
                        
                    end
                    
                end
            end
        end
        
    elseif  ( HEN == 4 )                                                    % Fully enriched element
                                                      
        if ( numel(PSI) == 0 )
            PN = [0 0 0 0]; 
        else
            PN = [ PSI{NN(1,4)}(N1)  PSI{NN(2,4)}(N2)  PSI{NN(3,4)}(N3)  PSI{NN(4,4)}(N4)];                 % Nodal crack level set values
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
            Nenr = [];
            
            Benru = [];
            Benrp = [];
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            index = 1;
            R_weight = 0;
            for iN = 1:4
                if NN(iN,2) ~= 0
                    icrack =  NN(iN,4);
                    psi1 = PSI{icrack}(N1);                                         % Psi level set value at node 1
                    psi2 = PSI{icrack}(N2);                                         % Psi level set value at node 2
                    psi3 = PSI{icrack}(N3);                                         % Psi level set value at node 3
                    psi4 = PSI{icrack}(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
    
                    Hgp = sign(psi);                                        % Heaviside value at current gauss point
                    Hi  = NN(iN,3);                                         % Nodal Heaviside value
                    H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value
                    
                    R_weight = R_weight + N(iN);
                    
                    Hgp = abs(psi);                                        % Heaviside value at current gauss point
                    Hi  = abs(PSI{icrack}(NN(iN,1)));                                         % Nodal Heaviside value
                    D   = (Hgp-Hi);                                           % Shifted distance value

                    Benru(:,(2*index-1):(2*index)) = [Nxy(1,iN)*H    0;
                                                     0              Nxy(2,iN)*H;
                                                     Nxy(2,iN)*H    Nxy(1,iN)*H];
                    Benrp(:,index) = [Nxy(1,iN)*D;Nxy(2,iN)*D];
                    Nenr(:,index) = N(1,iN)*D;
                    index = index+1;
                    if (i == 1)
                        localD(iLocD:(iLocD+1)) = [3*NN(iN,2)-2 3*NN(iN,2)-1];
                        iLocD = iLocD+2;
                        localP(iLocP) = 3*NN(iN,2);
                        iLocP = iLocP+1;
                    end
                end    
            end

            Bu = [Bu Benru];
            Bp = [Nxy Benrp.*R_weight];
            Nen = [N Nenr.*R_weight];
            DSPD=DISPDD(localD);
            deps = Bu*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            Strain = STATEV{iElem}{i}.strain;
            strainN = Strain + DEPS ;
            
            ELEPORP=DISPTD(localP);
            K = HydroConductivity_tan_update(PROP);
            velocityN = -K*Bp*ELEPORP;
            GPporepressure = Nen*ELEPORP;
            
            if PROP.nonlocal
                nonlocal_table = STATEV{iElem}{i}.nonlocalTable;
                [NLEquivStrain,scale] = computeNonlocalEquivalentStrain(nonlocal_table);
            else
                NLEquivStrain = computeEquivalentStrain(strainN);
            end
            
            Damage = STATEV{iElem}{i}.damage;
            Kappa = STATEV{iElem}{i}.kappa;
            [stressN, damageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,Kappa,Damage,GPporepressure);
                        
            % Update plastic variables once convergence is reached within one increment
            if UPDATE
                %%%%%%
                STATEV{iElem}{i}.sigma = stressN;
                STATEV{iElem}{i}.damage = damageN;
                STATEV{iElem}{i}.strain = strainN;
                STATEV{iElem}{i}.kappa = kappaN;
                STATEV{iElem}{i}.NLEquivStrain = NLEquivStrain;
                STATEV{iElem}{i}.fluidVelocity = velocityN;
                continue;
            end

            % Update residual Force by minus internal forces
            FORCE(localD) = FORCE(localD) - W*detJ*Bu'*[stressN(1:2,1); stressN(4,1)];

            M_coef=BiotModulus_tan_update(PROP);
            H_coef=HydroConductivity_tan_update(PROP);
            Q_coef=BiotCoefficient_tan_update(PROP);

            M = W*Nen'*M_coef*Nen*detJ;
            H = W*Bp'*H_coef*Bp*detJ;
            Q = W*Bu'*Q_coef'*Nen*detJ;
            QT = W*Nen'*Q_coef*Bu*detJ;
            
%             FORCE(localP) = FORCE(localP)+QT*( 1/(Theta*deltaT)*PREFTD(localD)+(1/Theta-1)*DISPTD_Rate(localD) )...
%                                           +M*( 1/(Theta*deltaT)*PREFTD(localP)+(1/Theta-1)*DISPTD_Rate(localP) )...
%                                           -1/(Theta*deltaT)*QT*DISPTD(localD) - (1/(Theta*deltaT)*M+H)*DISPTD(localP);
            FORCE(localP) = FORCE(localP) - QT*(DISPTD(localD)-PREFTD(localD)) - M*(DISPTD(localP)-PREFTD(localP))...
                            -deltaT*H*(Theta*DISPTD(localP)+(1-Theta)*PREFTD(localP)) ;
            % Tangent stiffness
            if LTAN
                C=Gauss_tan_update(PROP,damageN);                                     % Damage dependent stiffness matrix
                GKF(localD,localD) = GKF(localD,localD) + W*Bu'*C*Bu*detJ;                 % Sum up all Gauss point contribution                
%                 GKF(localP,localP) = GKF(localP,localP) + 1/(Theta*deltaT)*M + H;                
%                 GKF(localD,localP) = GKF(localD,localP) - Q;
%                 GKF(localP,localD) = GKF(localP,localD) + 1/(Theta*deltaT)*QT;
                GKF(localP,localP) = GKF(localP,localP) + M + H*(Theta*deltaT);
                GKF(localD,localP) = GKF(localD,localP) - Q;
                GKF(localP,localD) = GKF(localP,localD) + QT;
                
                if NLTAN
                    lcoeff = LocalCoefficient(PROP,kappaN,NLEquivStrain);
                    if ( lcoeff == 0 )
                        continue;
                    end
                    [dC_domega] = Gauss_tan_derivative(PROP);
                    
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

                        local_nl  = [N1_nl*3-2 N1_nl*3-1 N2_nl*3-2 N2_nl*3-1 N3_nl*3-2 N3_nl*3-1 N4_nl*3-2 N4_nl*3-1];             % Traditional index locations
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
                                icrack =  NN_nl(iN,4);
                                psi1 = PSI{icrack}(N1_nl);                                              % Psi level set value at node 1
                                psi2 = PSI{icrack}(N2_nl);                                              % Psi level set value at node 2
                                psi3 = PSI{icrack}(N3_nl);                                              % Psi level set value at node 3
                                psi4 = PSI{icrack}(N4_nl);                                              % Psi level set value at node 4
                                psi  = N_nl(1)*psi1+N_nl(2)*psi2+N_nl(3)*psi3+N_nl(4)*psi4;     % Psi level set value at current gauss point
    
                                Hgp = sign(psi);                                        % Heaviside value at current nonlocal gauss point
                                Hi  = NN_nl(iN,3);                                         % Nodal Heaviside value
                                H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

                                Ba = [Nxy_nl(1,iN)*H    0;
                                      0                 Nxy_nl(2,iN)*H;
                                      Nxy_nl(2,iN)*H    Nxy_nl(1,iN)*H];
                                  
                                Benr_nl(:,(2*index-1):(2*index)) = Ba;
                                index = index+1;
                                local_nl(iLoc_nl:(iLoc_nl+1)) = [3*NN_nl(iN,2)-2 3*NN_nl(iN,2)-1];
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
                        if (equivalentEPS_nl ~= 0 && lcoeff ~= 0)
                            [esp_pri_dir,esp_pri_val]= eigs([strain_nl(1,1) strain_nl(4,1); strain_nl(4,1) strain_nl(2,1)]);
                            angle_eps_pri=[0;0];
                            if (esp_pri_val(1,1) > 0)
                                angle_eps_pri(1) = esp_pri_val(1,1);
                            end
                            if (esp_pri_val(2,2) > 0)
                                angle_eps_pri(2) = esp_pri_val(2,2);
                            end
                            dequieps_deps = [esp_pri_dir(1,1)*angle_eps_pri(1)   esp_pri_dir(2,2)*angle_eps_pri(2)  esp_pri_dir(1,2)*angle_eps_pri(2)+esp_pri_dir(2,1)*angle_eps_pri(1)]./equivalentEPS_nl;
                            nl_contribution = (weight_nl*volume_nl/scale).*(dequieps_deps*B_nl);
                            CBu_nl = [strainN(1:2,1); 2*strainN(4,1)]*nl_contribution;
                            GKF(localD,local_nl) = GKF(localD,local_nl) + lcoeff*detJ*W*Bu'*dC_domega*CBu_nl;
                        end
                    end
                end
            end
        end
        
        for icrack = 1:nCrack
            I_frac = find( connec_frac{icrack}(:,1) == iElem );
            if ~isempty(I_frac)
                break;
            end
        end
        
        if isempty(I_frac)
            continue;
        end
        
        N1_frac  = connec_frac{icrack}(I_frac,2);                                              % Node 1 for current fracture segment
        N2_frac  = connec_frac{icrack}(I_frac,3);                                              % Node 2 for current fracture segment
        
        X1_frac = xyz_frac{icrack}(N1_frac,4);
        X2_frac = xyz_frac{icrack}(N2_frac,4);
        
        tan_frac = xyz_frac{icrack}(N2_frac,2:3) - xyz_frac{icrack}(N1_frac,2:3);
        tan_frac = tan_frac./sqrt(sum(tan_frac.^2));
        norm_frac = [-tan_frac(2) tan_frac(1)];
        
        lambda = [tan_frac; norm_frac];
        
%         local_frac  = [xyz_frac(N1_frac,1) xyz_frac(N2_frac,1)];            % Traditional index locations
        localD_enr  = localD(1,9:16);
        localD_une  = localD(1,1:8);
        
        Ngp = size(statev_frac{icrack}{I_frac},2);
        
        for ig = 1:Ngp

%             gp_frac = statev_frac{I_frac}{ig}.natural_coodinates;            % Gauss points
            gp_ele = statev_frac{icrack}{I_frac}{ig}.natural_coordinate_ele;            % Gauss points
            W = statev_frac{icrack}{I_frac}{ig}.gauss_weight;                        % Gauss weights
%             czmlength = statev_frac{I_frac}{ig}.seglength;  
            histo_jump = statev_frac{icrack}{I_frac}{ig}.jump;
            frac_prop = statev_frac{icrack}{I_frac}{ig};
            
            xi = gp_ele(1,1); eta = gp_ele(2,1);                                % Gauss points
            
            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                     % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];
                  
            [Nxy, ~] = Shape_Function(xi, eta, xyz);
                  
%             N_frac = 1/2*[(1-gp_frac)   (1+gp_frac)];
            detJ_frac = det([-0.5 0.5]*[X1_frac; X2_frac]);
%             Nxy_frac = 1/detJ_frac*[-0.5 0.5];
 
            Nenr = [];
            Benrp = [];

            index = 1;
            R_weight = 0;
            for iN = 1:4
%                 Hi  = NN(iN,3);                                         % Nodal Heaviside value
%                 H   = (0-Hi)/2;

                R_weight = R_weight + N(iN);
                    
                Hi  = abs(PSI{icrack}(NN(iN,1)));                                         % Nodal Heaviside value
                D   = (0-Hi);                                           % Shifted distance value

%                 Benru(:,(2*index-1):(2*index)) = [Nxy(1,iN)*H    0;
%                                                   0              Nxy(2,iN)*H;
%                                                   Nxy(2,iN)*H    Nxy(1,iN)*H];
                Benrp(:,index) = [Nxy(1,iN)*D;Nxy(2,iN)*D];
                Nenr(:,index) = N(1,iN)*D;
                index = index+1;
            end
            
            Bp = [Nxy Benrp.*R_weight];
            Nen = [N Nenr.*R_weight];
            NT = [N(1)   0     N(2)  0     N(3)  0     N(4)  0 ;
                  0      N(1)  0     N(2)  0     N(3)  0     N(4);];
            
            jump = lambda*NT*DISPTD(localD_enr);
            fracwidth = jump(2);
            
            if (histo_jump(1)==1e6)
                if fracwidth<1
                    fracwidth=1;
                end
                Traction=[0;0];
            else
                Traction =  cohesiveLaw(jump, histo_jump, frac_prop);
            end
                
            new_jump=histo_jump;
            if jump(1) >= histo_jump(1)
                new_jump(1) = jump(1);
            end
            if jump(2) >= histo_jump(2)
                new_jump(2) = jump(2);
            end
            
            if UPDATE
                %%%%%%
                statev_frac{icrack}{I_frac}{ig}.width = fracwidth;
                statev_frac{icrack}{I_frac}{ig}.jump = new_jump;
                continue;
            end
            
            %% Update residual Force by minus internal forces    
            
            L_up =  W*detJ_frac*NT'*norm_frac'*Nen;
            L_upT = W*detJ_frac*Nen'*norm_frac*NT;
          
            M_coef=FracBiotModulus_tan_update(PROP,fracwidth);
            K_coef=FracHydroConductivity_tan_update(PROP,fracwidth);      
            M_frac = W*Nen'*M_coef*Nen*detJ_frac;
            H_frac = W*Bp'*tan_frac'*K_coef*tan_frac*Bp*detJ_frac;   
            
            FORCE(localD_enr) = FORCE(localD_enr) - W*detJ_frac*NT'*lambda'*Traction... 
                                                  + W*detJ_frac*NT'*norm_frac'*Nen*DISPTD(localP);
%             FORCE(localD_une) = FORCE(localD_une) - W*detJ_frac*NT'*norm_frac'*Nen*DISPTD(localP);
                                              
            FORCE(localP) = FORCE(localP) - L_upT*(DISPTD(localD_enr)-PREFTD(localD_enr))...
                                          - M_frac*(DISPTD(localP)-PREFTD(localP))...
                                          - deltaT*H_frac*(Theta*DISPTD(localP)+(1-Theta)*PREFTD(localP));
            
%             FORCE(localP) = FORCE(localP) - W*detJ_frac*Nen'*M_coef*Nen*(1/(Theta*deltaT)*(DISPTD(localP)-PREFTD(localP))-(1/Theta-1)*DISPTD_Rate(localP))...
%                                           - W*detJ_frac*Nen'*(norm_frac*NT*( 1/(Theta*deltaT)*(DISPTD(localD_enr)-PREFTD(localD_enr))+(1/Theta-1)*DISPTD_Rate(localD_enr)))...
%                                           - W*detJ_frac*Bp'*tan_frac'*K_coef*tan_frac*Bp*DISPTD(localP);
%                                           - W*detJ_frac*Nen'*M_coef*(1/(Theta*deltaT)*Nen*DISPTD(localP))...
%                                           - W*detJ_frac*Nen'*(norm_frac*( 1/(Theta*deltaT)*NT*DISPTD(localD_enr)))...
                                          
            if LTAN
                if (histo_jump(1)==1e6)
                    M_czm=[0 0;0 0;];
                else
                    M_czm =  cohesiveStiffness(jump, histo_jump, frac_prop);
                end
                K_CZM = W*NT'*lambda'*M_czm*lambda*NT*detJ_frac;  %% cohesive zone section
                
                GKF(localD_enr,localD_enr) = GKF(localD_enr,localD_enr) + K_CZM;  
%                 GKF(localD_enr,localP) = GKF(localD_enr,localP) + L_up;
%                 GKF(localP,localD_enr) = GKF(localP,localD_enr) - L_up';
%                 GKF(localP,localP) = GKF(localP,localP) - (1/(Theta*deltaT)*M_frac + H_frac);   
                GKF(localD_enr,localP) = GKF(localD_enr,localP) - L_up;
                GKF(localP,localD_enr) = GKF(localP,localD_enr) + L_upT;
                GKF(localP,localP) = GKF(localP,localP) + M_frac + (Theta*deltaT)*H_frac; 
                
            end
        end

    end 

end

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
    

    function [stressN, DamageN, kappaN]=Gauss_sig_update(PROP,NLEquivStrain,strainN,kappa,Damage,GPporepressure)
        % Inputs:
        % PROP = [E0, nu0, epsilon_0, epsilon_f];
        % D = 4*4 elastic stiffness matrix
        % stress = [s11, s22, s12,   s33];
        % strain = [e11, e22, 2*e12, e33];
        % Plane strain problem e33=0
        %%
        
            if ( NLEquivStrain > kappa )
                kappa = NLEquivStrain;
                Damage = 1-PROP.eps_cr/kappa*exp(-PROP.B*(kappa - PROP.eps_cr));
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
        
            alpha = BiotCoefficient_tan_update(PROP);
        
            stressN = MATC*[strainN(1:3,1); 2*strainN(4,1)] - [alpha(1);alpha(2);alpha(1);alpha(3)]* GPporepressure;    %updated stress
        
    end

    function [MATC]=Gauss_tan_update(PROP,Damage)
    % Inputs:
    % Plane strain problem e33=0
        MATC=zeros(3,3);
        D= (1-Damage)*PROP.E/(1+PROP.nu)/(1-2*PROP.nu);
        MATC(1,1)=(1-PROP.nu)*D;
        MATC(1,2)=PROP.nu*D;
        MATC(1,3)=0;
        MATC(2,1)=PROP.nu*D;
        MATC(2,2)=(1-PROP.nu)*D;
        MATC(2,3)=0;
        MATC(3,1)=0;
        MATC(3,2)=0;
        MATC(3,3)=(1-2*PROP.nu)/2*D;
    end

    function [dC_domega]=Gauss_tan_derivative(PROP)
        % Inputs: derivative stiffness matrix with respective to damage scalar
            dC_domega = zeros(3,3);
            D= -PROP.E/(1+PROP.nu)/(1-2*PROP.nu);
            dC_domega(1,1)=(1-PROP.nu)*D;
            dC_domega(1,2)=PROP.nu*D;
            dC_domega(1,3)=0;
            dC_domega(2,1)=PROP.nu*D;
            dC_domega(2,2)=(1-PROP.nu)*D;
            dC_domega(2,3)=0;
            dC_domega(3,1)=0;
            dC_domega(3,2)=0;
            dC_domega(3,3)=(1-2*PROP.nu)/2*D;
    end

function [traction] =  cohesiveLaw(jump, histo_jump, frac_prop)

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
        traction(2)=frac_prop.PenaltyStiffness*deln;
    elseif( deln >= 0 && deln <= frac_prop.deltaN  && delt <= frac_prop.deltaT_conj )
        if (deln >= deln_max )
            traction(2)=(frac_prop.GammaN/frac_prop.deltaN)*( (frac_prop.m*(1-deln/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m)*...
                                (frac_prop.GammaT*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n+frac_prop.dGtn);
        else
            traction(2)=(frac_prop.GammaN/frac_prop.deltaN)*( (frac_prop.m*(1-deln_max/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln_max/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^frac_prop.m)*...
                                (frac_prop.GammaT*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n+frac_prop.dGtn)*deln/deln_max;
        end
    end
    if ( deln >= 0 && deln <= frac_prop.deltaN_conj && delt <= frac_prop.deltaT )
        if ( delt >= delt_max )
            traction(1)=(frac_prop.GammaT/frac_prop.deltaT)*( (frac_prop.n*(1-delt/frac_prop.deltaT)^frac_prop.beta)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^(frac_prop.n-1)-...
                                frac_prop.beta*(1-delt/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n)*...
                                (frac_prop.GammaN*(1-deln/frac_prop.deltaN)^frac_prop.alpha*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m+frac_prop.dGtn)*sign_dt;  
        else
            traction(1)=(frac_prop.GammaT/frac_prop.deltaT)*( (frac_prop.n*(1-delt_max/frac_prop.deltaT)^frac_prop.beta)*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^(frac_prop.n-1)-...
                                frac_prop.beta*(1-delt_max/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^frac_prop.n)*...
                                (frac_prop.GammaN*(1-deln/frac_prop.deltaN)^frac_prop.alpha*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m+frac_prop.dGtn)*delt*sign_dt/delt_max;
        end
    end
end


function [stiffness] =  cohesiveStiffness(jump, histo_jump, frac_prop)

%     global frac_prop
    
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
        stiffness(2,2)=frac_prop.PenaltyStiffness;
        stiffness(2,1)=0;
    elseif( deln >= 0 && deln <= frac_prop.deltaN  && delt <= frac_prop.deltaN_conj )
        if (deln >= deln_max )
            stiffness(2,2)=(frac_prop.GammaN/frac_prop.deltaN^2)*( ((frac_prop.m^2-frac_prop.m)*(1-deln/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^(frac_prop.m-2)+...
                                ((frac_prop.alpha^2-frac_prop.alpha)*(1-deln/frac_prop.deltaN)^(frac_prop.alpha-2)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m)-...
                                2*frac_prop.alpha*frac_prop.m*(1-deln/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^(frac_prop.m-1))*...
                                (frac_prop.GammaT*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n+frac_prop.dGtn);
            stiffness(2,1)=(frac_prop.GammaN*frac_prop.GammaT/frac_prop.deltaN/frac_prop.deltaT)*( ( frac_prop.m*(1-deln/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m)*...
                                ( frac_prop.n*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^(frac_prop.n-1)-...
                                  frac_prop.beta*(1-delt/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n )*sign_dt;
        else
            stiffness(2,2)=(frac_prop.GammaN/frac_prop.deltaN)*( (frac_prop.m*(1-deln_max/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln_max/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^frac_prop.m)*...
                                (frac_prop.GammaT*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n+frac_prop.dGtn)/deln_max;
            stiffness(2,1)=(frac_prop.GammaN*frac_prop.GammaT/frac_prop.deltaN/frac_prop.deltaT)*( ( frac_prop.m*(1-deln_max/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln_max/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln_max/frac_prop.deltaN)^frac_prop.m)*...
                                ( frac_prop.n*(1-delt/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^(frac_prop.n-1)-...
                                  frac_prop.beta*(1-delt/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n )*sign_dt*deln/deln_max;
        end
    end
    if ( deln >= 0 && deln <= frac_prop.deltaN_conj && delt <= frac_prop.deltaT )
        if ( delt >= delt_max )
            stiffness(1,1)=(frac_prop.GammaT/frac_prop.deltaT^2)*( ((frac_prop.n^2-frac_prop.n)*(1-delt/frac_prop.deltaT)^frac_prop.beta)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^(frac_prop.n-2)+...
                                ((frac_prop.beta^2-frac_prop.beta)*(1-delt/frac_prop.deltaT)^(frac_prop.beta-2)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^frac_prop.n)-...
                                2*frac_prop.beta*frac_prop.n*(1-delt/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt/frac_prop.deltaT)^(frac_prop.n-1))*...
                                (frac_prop.GammaN*(1-deln/frac_prop.deltaN)^frac_prop.alpha*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.n+frac_prop.dGnt);
            stiffness(1,2)=stiffness(2,1);
        else
            stiffness(1,1)=(frac_prop.GammaT/frac_prop.deltaT)*( (frac_prop.n*(1-delt_max/frac_prop.deltaT)^frac_prop.beta)*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^(frac_prop.n-1)-...
                                frac_prop.beta*(1-delt_max/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^frac_prop.n)*...
                                (frac_prop.GammaN*(1-deln/frac_prop.deltaN)^frac_prop.alpha*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m+frac_prop.dGtn)/delt_max;
            stiffness(1,2) = (frac_prop.GammaN*frac_prop.GammaT/frac_prop.deltaN/frac_prop.deltaT)*( ( frac_prop.m*(1-deln/frac_prop.deltaN)^frac_prop.alpha)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^(frac_prop.m-1)-...
                                frac_prop.alpha*(1-deln/frac_prop.deltaN)^(frac_prop.alpha-1)*(frac_prop.m/frac_prop.alpha+deln/frac_prop.deltaN)^frac_prop.m)*...
                                ( frac_prop.n*(1-delt_max/frac_prop.deltaT)^frac_prop.beta*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^(frac_prop.n-1)-...
                                  frac_prop.beta*(1-delt_max/frac_prop.deltaT)^(frac_prop.beta-1)*(frac_prop.n/frac_prop.beta+delt_max/frac_prop.deltaT)^frac_prop.n )*sign_dt*delt/delt_max;    
        end
    end

end

function [lcoeff] = LocalCoefficient(PROP,kappa,NLEquivStrain)
    % the local derivative of damage with respect to kappa
    % demage eovlution function is defined as
    % Damage = 1-PROP.eps_cr/kappa*exp(-PROP.B*(kappa - PROP.eps_cr));
        if (NLEquivStrain == kappa) 
            lcoeff = (1/kappa + PROP.B)*(PROP.eps_cr/kappa)*exp(-PROP.B*(kappa-PROP.eps_cr));
        else
            lcoeff =  0;
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


function [M] = FracBiotModulus_tan_update(PROP,width)
    M = width/PROP.Kf;
end

function [K] = FracHydroConductivity_tan_update(PROP,width)
    K = width^3/12/PROP.viscosity;
end

function [C] = FracLeakoffCoef_tan_update(PROP)
    C = PROP.leakoff;
end

function [M]=BiotModulus_tan_update(PROP)

	Kf = PROP.Kf;                   %mm^2
    phi = PROP.porosity;            %mm^2
    alpha = mean(PROP.BiotAlpha);
    M = phi/Kf +(alpha-phi)/PROP.Ks;

end

function [Q]=BiotCoefficient_tan_update(PROP)

    Q = [PROP.BiotAlpha 0];
    
end

function [K] = HydroConductivity_tan_update(PROP) 
% Inputs:
% Plane strain problem e33=0
%%
    kappa11 = PROP.kappa11;           %mm^2
    kappa22 = PROP.kappa22;            %mm^2
    viscosity = PROP.viscosity;     %MPa*s            0.0010005 Pa*s

    K=zeros(2,2);

    K(1,1) = kappa11/viscosity;
    K(2,2) = kappa22/viscosity;

end

function [EquivStrain] = computeEquivalentStrain( strain)
    %%
        [~,strain_principal]= eigs([strain(1,1) strain(4,1); strain(4,1) strain(2,1)]);
        EquivStrain = sqrt((strain_principal(1,1)+abs(strain_principal(1,1))/2)^2+...
                           (strain_principal(2,2)+abs(strain_principal(2,2))/2)^2);   
end