% Written By: Wencheng Jin, Idaho National Laboratory (2020)
% Website: https://sites.google.com/view/wenchengjin/software
% Email: wencheng.jin@inl.gov

function updateDomainBeforeNonlocAverage(PROP)
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

global CONNEC NODES PSI XYZ STATEV DISPDD

% Initialize the FE stiffness matrix
% if iter > 1, globalDOF = globalDOF+16*iter; end                             % Initialize extra space for growing
% GKF = sparse(globalDOF,globalDOF);                                      % Define the global K

thickness = PROP.plane_thickness;

for iElem = 1:size(CONNEC,1)
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
    HEN = nnz(NN(:,2));
    SIGN = abs(sum(NN(:,3)));% Number of nodes with Heaviside enrichment
    local  = [N1*3-2 N1*3-1 N2*3-2 N2*3-1 N3*3-2 N3*3-1 N4*3-2 N4*3-1];             % Traditional index locations

    iLoc   = 9; 
    
    if (HEN == 0)                                                           % Unenriched nodes
        % Traditional element
        X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2); % Nodal x-coordinates
        Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3); % Nodal y-coordinates     
        xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix   
        
        [gp,gw] = gauss(2,'QUAD');
        for i = 1:length(gp)
           xi = gp(i,1); eta = gp(i,2);                                % Gauss points
           [Nxy, detJ] = Shape_Function(xi, eta, xyz);                % Derivative of shape functions with respect to x/y;  Determinant of the Jacobian
           
           Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];
           
           DSPD=DISPDD(local);
           deps = Bu*DSPD;
           DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];
            
           Strain = STATEV{iElem}{i}.strain;
            
           Strain=Strain+DEPS;
           EquivStrain = computeEquivalentStrain(Strain); 
            
           STATEV{iElem}{i}.EquivStrain = EquivStrain;
           STATEV{iElem}{i}.volume = detJ*gw(i,1)*thickness;
        end         
        
    elseif HEN > 0                                                                 % Enriched element
        
        X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
        Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
        xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix

        if HEN == 4 && SIGN ~= 4                                                    % Fully enriched element
            if numel(PSI) == 0, PN = [0 0 0 0]; else
                  icrack = NN(1,4);
                  PN = [ PSI{icrack}(N1)  PSI{icrack}(N2)  PSI{icrack}(N3)  PSI{icrack}(N4)];               % Nodal crack level set values
            end
            [gp,gw,J] = subDomain(3,PN,xyz);                         % Full Heaviside enrichment
        else                                                                % Partially enriched element
            [gp,gw] = gauss(2,'QUAD');
            J = [];
        end

        for i = 1:length(gp)
            xi = gp(i,1); eta = gp(i,2);                                    % Gauss points
            if isempty(J) == 0
                Ji   = [J(i,1) J(i,2);J(i,3) J(i,4)];                       % Jacobian of subdomain
                detJ = det(Ji);                                             % Determinant of the Jacobian
                [Nxy, ~] = Shape_Function(xi, eta, xyz);
            else
                [Nxy, detJ] = Shape_Function(xi, eta, xyz);                 % Determinant of the Jacobian and the derivative of shape functions with respect to x/y                                 
            end
            
            N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...                     % Shape functions
                      (1+xi)*(1+eta);(1-xi)*(1+eta)];

            Benr = [];
            Bu = [Nxy(1,1)   0          Nxy(1,2)   0          Nxy(1,3)   0          Nxy(1,4)   0;...
                  0          Nxy(2,1)   0          Nxy(2,2)   0          Nxy(2,3)   0          Nxy(2,4);...
                  Nxy(2,1)   Nxy(1,1)   Nxy(2,2)   Nxy(1,2)   Nxy(2,3)   Nxy(1,3)   Nxy(2,4)   Nxy(1,4)];

            index = 1;
            for iN = 1:4
                if NN(iN,2) ~= 0
                    icrack = NN(iN,4);
                    psi1 = PSI{icrack}(N1);                                         % Psi level set value at node 1
                    psi2 = PSI{icrack}(N2);                                         % Psi level set value at node 2
                    psi3 = PSI{icrack}(N3);                                         % Psi level set value at node 3
                    psi4 = PSI{icrack}(N4);                                         % Psi level set value at node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
    
                    Hgp = sign(psi);                                        % Heaviside value at current gauss point
                    Hi  = NN(iN,3);                                         % Nodal Heaviside value
                    H   = (Hgp-Hi)/2;                                           % Shifted Heaviside value

                    Ba = [Nxy(1,iN)*H    0;
                          0              Nxy(2,iN)*H;
                          Nxy(2,iN)*H    Nxy(1,iN)*H];
                    Benr(:,(2*index-1):(2*index)) = Ba;
                    index = index+1;
                    if (i == 1)
                        local(iLoc:(iLoc+1)) = [3*NN(iN,2)-2 3*NN(iN,2)-1];
                        iLoc = iLoc+2;
                    end
                end    
            end

            B = [Bu Benr];
 
            DSPD=DISPDD(local);
            deps = B*DSPD;
            DEPS =[deps(1,1); deps(2,1); 0; 0.5*deps(3,1);];

            Strain = STATEV{iElem}{i}.strain;
            
            Strain=Strain+DEPS;
            EquivStrain = computeEquivalentStrain(Strain); 
            
            STATEV{iElem}{i}.EquivStrain = EquivStrain;
            STATEV{iElem}{i}.volume = detJ*gw(i,1)*thickness;
        end

    end 

end

% GKF = GKF + sparse(allRowT,allColT,allValT,globalDOF,globalDOF);

end

function [EquivStrain] = computeEquivalentStrain( strain)
    %%
        [~,strain_principal]= eigs([strain(1,1) strain(4,1); strain(4,1) strain(2,1)]);
        EquivStrain = sqrt((strain_principal(1,1)+abs(strain_principal(1,1))/2)^2+...
                           (strain_principal(2,2)+abs(strain_principal(2,2))/2)^2);   
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

