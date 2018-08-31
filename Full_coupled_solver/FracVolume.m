% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function FracVolume()
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

global CONNEC NODES DISPTD connec_frac xyz_frac statev_frac Fracvolume DispJump

enrElemwithTip = enrElem;
m=size(enrElemwithTip,1);
EnrichElems = enrElemwithTip(1:m,1);

DispJump = zeros(m+1,1);
% check = 0;

% WIDTH = zeros(size(PRESF));
% seglenth = zeros(m-1,1);

for ie = 1:m
    iElem = EnrichElems(ie);
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
    HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    
    localD  = [N1*3-2 N1*3-1 N2*3-2 N2*3-1 N3*3-2 N3*3-1 N4*3-2 N4*3-1];    % Traditional index locations for displacement
    iLocD   = 9;                                                    % Fully enriched element
                                                      
        I_frac = find( connec_frac(:,1) == iElem );
        
        N1_frac  = connec_frac(I_frac,2);                                              % Node 1 for current fracture segment
        N2_frac  = connec_frac(I_frac,3);                                              % Node 2 for current fracture segment
        
        X1_frac = xyz_frac(N1_frac,4);
        X2_frac = xyz_frac(N2_frac,4);
        
        
        local_frac  = [xyz_frac(N1_frac,1) xyz_frac(N2_frac,1)];
        
        tan_frac = xyz_frac(N2_frac,2:3) - xyz_frac(N1_frac,2:3);
        tan_frac = tan_frac./sqrt(sum(tan_frac.^2));
        norm_frac = [-tan_frac(2) tan_frac(1)];
                
        Ngp = size(statev_frac{I_frac},2);
        
%         meanwidth = 0;
        
        for ig = 1:Ngp

            gp_frac = statev_frac{I_frac}{ig}.natural_coodinates;            % Gauss points
            gp_ele = statev_frac{I_frac}{ig}.natural_coordinate_ele;            % Gauss points
            W = statev_frac{I_frac}{ig}.gauss_weight;                        % Gauss weights
%             czmlength = statev_frac{I_frac}{ig}.seglength;
            
            xi = gp_ele(1,1); eta = gp_ele(2,1);                                % Gauss points
            
            N  = 1/4*[(1-xi)*(1-eta) (1+xi)*(1-eta) ...                     % Shape functions
                      (1+xi)*(1+eta) (1-xi)*(1+eta)];
                  
            N_frac = 1/2*[(1-gp_frac)   (1+gp_frac)];
            detJ_frac = det([-0.5 0.5]*[X1_frac; X2_frac]);
            
            for iN = 1:4
                localD(iLocD:(iLocD+1)) = [3*NN(iN,2)-2 3*NN(iN,2)-1];
                iLocD = iLocD+2;
            end
            
            localD_enr  = localD(1,9:16);
            
            NT = [N(1)  0     N(2)  0     N(3)  0     N(4)    0;
                  0     N(1)  0     N(2)  0     N(3)  0       N(4);];
            
            fracwidth = norm_frac*NT*DISPTD(localD_enr);
%             meanwidth =  meanwidth + fracwidth;
%             delta_fracwidth = norm_frac*NT*DISPDD(localD_enr);
            
            DispJump(local_frac) = DispJump(local_frac) + W*detJ_frac*N_frac'*fracwidth;
                                       
        end
        
%         check = check + czmlength*meanwidth/2;
end
Fracvolume = sum(DispJump);
end



