% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function Postprocessor(STEP,EXTDISP,file_name)
% This function output the all the primary and secondary results of the
% simulation into vtk formate for post processing

global  PHI PSI XYZ CONNEC NODES STATEV DISPTD connec_frac xyz_frac statev_frac FORCE

%Recover the state varibals at the nodes from Gauss Points

psi = full(PSI);
phi = full(PHI);

DISPTD_ = full(DISPTD);

npoints=size(XYZ,1);
ncells=size(CONNEC,1);
if isempty(psi)
    psi = zeros(npoints,1);
    phi = zeros(npoints,1);
end
stress = zeros(npoints,4);
damage = zeros(npoints,1);
strain = zeros(npoints,4);
NLEquivStrain = zeros(npoints,1);

EnrichElems = enrElem();
UnenrichedElem = setdiff([1:ncells]', EnrichElems);
NELE=size(EnrichElems,1);

%
%Extrapolation for regular elememt nodes
% EnrichElems = [];
% for ielem=1:ncells
%     if size(STATEV{ielem},2) ~= 4
%         EnrichElems = [EnrichElems; ielem;]; 
%         continue;
%     end
%     Nodes = CONNEC(ielem,2:5)';
%     A = [1+sqrt(3)/2   -0.5          1-sqrt(3)/2    -0.5;
%          -0.5          1+sqrt(3)/2   -0.5           1-sqrt(3)/2;
%          1-sqrt(3)/2   -0.5          1+sqrt(3)/2    -0.5;
%          -0.5          1-sqrt(3)/2   -0.5           1+sqrt(3)/2;];
%     Sigma = zeros(4,4);
%     for j=1:4
%         Sigma(:,j)= A*[STATEV{ielem}{1}.sigma(j,1); STATEV{ielem}{2}.sigma(j,1); STATEV{ielem}{3}.sigma(j,1); STATEV{ielem}{4}.sigma(j,1)];
%     end
%     Epsilon = zeros(4,4);
%     for j=1:4
%         Epsilon(:,j)= A*[STATEV{ielem}{1}.strain(j,1); STATEV{ielem}{2}.strain(j,1); STATEV{ielem}{3}.strain(j,1); STATEV{ielem}{4}.strain(j,1)];
%     end
%     Omega = zeros(4,2);
%     for j=1:2
%         Omega(:,j)= A*[STATEV{ielem}{1}.damage(j,1); STATEV{ielem}{2}.damage(j,1); STATEV{ielem}{3}.damage(j,1); STATEV{ielem}{4}.damage(j,1)];
%     end
%     Equ_epsilon = zeros(4,2);
%     for j=1:2
%         Equ_epsilon(:,j)= A*[STATEV{ielem}{1}.NLEquivStrain(j,1); STATEV{ielem}{2}.NLEquivStrain(j,1); STATEV{ielem}{3}.NLEquivStrain(j,1); STATEV{ielem}{4}.NLEquivStrain(j,1)];
%     end
%     stress(Nodes,:) = stress(Nodes,:) + Sigma;
%     damage(Nodes,:) = damage(Nodes,:) + Omega;
%     strain(Nodes,:) = strain(Nodes,:) + Epsilon;
%     NLEquivStrain(Nodes,:) = NLEquivStrain(Nodes,:) + Equ_epsilon;
% end

% UnenrichedElem = setdiff([1:ncells]', EnrichElems);

% allnodes = sort([CONNEC(UnenrichedElem,2); CONNEC(UnenrichedElem,3); CONNEC(UnenrichedElem,4); CONNEC(UnenrichedElem,5);]);
% [a,~]=hist(allnodes,unique(allnodes));
% for i=1:4
%     stress(:,i) = stress(:,i)./a';
% end
% for i=1:4
%     strain(:,i) = strain(:,i)./a';
% end
% for i=1:2
%     damage(:,i) = damage(:,i)./a';
% end
% for i=1:2
%     NLEquivStrain(:,i) = NLEquivStrain(:,i)./a';
% end
%%
% 
% %%%%%Patch Recovery methods for the nodes with enriched degree of freedom
% I=find(NODES(:,2)~=0);
% m=size(I,1);
for i=1:npoints
%     NodeNumber = NODES(I(i),1);
    NodeNumber = i;
    I1 = find(CONNEC(:,2) == NodeNumber);
    I2 = find(CONNEC(:,3) == NodeNumber);
    I3 = find(CONNEC(:,4) == NodeNumber);
    I4 = find(CONNEC(:,5) == NodeNumber);
    hElem    = unique([I1; I2; I3; I4;]);                                       % Candidate elements connected to the nodes
    A=zeros(3,3);
    b_stress1=zeros(3,1);
    b_stress2=zeros(3,1);
    b_stress3=zeros(3,1);
    b_stress4=zeros(3,1);
    b_damage1 = zeros(3,1);
    b_NLEquivStrain1= zeros(3,1);
    b_strain1= zeros(3,1);
    b_strain2= zeros(3,1);
    b_strain4= zeros(3,1);
    
    [~, n] = size(hElem);
    for j=1:n
        m=size(STATEV{hElem(j)},2);
        for k=1:m
            P=[1 STATEV{hElem(j)}{k}.true_coodinates(1) STATEV{hElem(j)}{k}.true_coodinates(2)];
            A = A + P'*P;
            b_stress1 = b_stress1+P'.*STATEV{hElem(j)}{k}.sigma(1,1);
            b_stress2 = b_stress2+P'.*STATEV{hElem(j)}{k}.sigma(2,1);
            b_stress3 = b_stress3+P'.*STATEV{hElem(j)}{k}.sigma(3,1);
            b_stress4 = b_stress4+P'.*STATEV{hElem(j)}{k}.sigma(4,1);
            b_damage1 = b_damage1+P'.*STATEV{hElem(j)}{k}.damage(1,1);
            b_NLEquivStrain1 = b_NLEquivStrain1+P'.*STATEV{hElem(j)}{k}.NLEquivStrain(1,1);
            b_strain1 = b_strain1+P'.*STATEV{hElem(j)}{k}.strain(1,1);
            b_strain2 = b_strain2+P'.*STATEV{hElem(j)}{k}.strain(2,1);
            b_strain4 = b_strain4+P'.*STATEV{hElem(j)}{k}.strain(4,1);
        end
    end
    a_stress1 = A\b_stress1;
    a_stress2 = A\b_stress2;
    a_stress3 = A\b_stress3;
    a_stress4 = A\b_stress4;
    a_damage1 = A\b_damage1;
    a_NLEquivStrain1 = A\b_NLEquivStrain1;
    a_strain1 = A\b_strain1;
    a_strain2 = A\b_strain2;
    a_strain4 = A\b_strain4;
    x = XYZ(i,2);
    y = XYZ(i,3);
    P = [1 x y];
    stress(NodeNumber,1)=P*a_stress1;
    stress(NodeNumber,2)=P*a_stress2;
    stress(NodeNumber,3)=P*a_stress3;
    stress(NodeNumber,4)=P*a_stress4;
    damage(NodeNumber,1)=P*a_damage1; 
    strain(NodeNumber,1)=P*a_strain1;
    strain(NodeNumber,2)=P*a_strain2;
    strain(NodeNumber,4)=P*a_strain4;
    NLEquivStrain(NodeNumber,1)=P*a_NLEquivStrain1;
end
% %
% 
% 
filename = strcat(file_name,num2str(STEP),'.vtu');
fid = fopen(filename, 'w');
% fprintf(fid, '<!-- TimeStep %8.6e -->\n',TIME);
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '<UnstructuredGrid>\n');
fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', npoints,ncells-NELE);
fprintf(fid, '<Points>\n');
fprintf(fid, '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, ' %8.6e %8.6e %8.6e', XYZ(i,2), XYZ(i,3), 0);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '</Points>\n');
fprintf(fid, '<Cells>\n');
fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
for i = 1:size(UnenrichedElem,1)
    iEle =  UnenrichedElem(i,1);
    fprintf(fid, ' %d %d %d %d', CONNEC(iEle,2)-1, CONNEC(iEle,3)-1, CONNEC(iEle,4)-1,CONNEC(iEle,5)-1);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
for i = 1:(ncells-NELE)
    fprintf(fid, ' %d', i*4);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
for i = 1:(ncells-NELE)
    fprintf(fid, ' %d', 9);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, ' </Cells>\n');
% fprintf(fid, '<PointData Scalars="LS_psi LS_phi " Vectors="DisplacementVector " Tensors="StressTensor DamageTensor " >\n');
% fprintf(fid, '<PointData Scalars="LS_psi LS_phi" Vectors="DisplacementVector StressTensor DamageTensor EquivalentStrain">\n');
fprintf(fid, '<PointData Scalars="LS_psi LS_phi" Vectors="Displacement Stress Strain Damage NLEquivalentStrain">\n');
fprintf(fid, '<DataArray type="Float64" Name="Displacement" NumberOfComponents="3" format="ascii">\n');
for i = 1:npoints
%     if (NODES(i,2)~=0)
%         fprintf(fid, '  %8.6e   %8.6e   %8.6e', DISPTD_(2*i-1)+DISPTD_(2*NODES(i,2)-1),DISPTD_(2*i)+DISPTD_(2*NODES(i,2)),0);
%     else
        fprintf(fid, '  %8.6e   %8.6e   %8.6e', DISPTD_(2*i-1),DISPTD_(2*i-0),0);
%     end
end 
fprintf(fid, ' </DataArray>\n');
% fprintf(fid, '<DataArray type="Float64" Name="PorePressure" NumberOfComponents="1" format="ascii">\n');
% for i = 1:npoints
%         fprintf(fid, '  %8.6e', PRESTD_(i));
% end 
% fprintf(fid, ' </DataArray>\n');
% fprintf(fid, '<DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
% for i = 1:npoints
%         fprintf(fid, '  %8.6e   %8.6e   %8.6e', Velocity(i,1), Velocity(i,2), 0);
% end 
% fprintf(fid, ' </DataArray>\n');


fprintf(fid, '<DataArray type="Float64" Name="Stress" NumberOfComponents="4" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e   %8.6e   %8.6e   %8.6e', stress(i,1), stress(i,2), stress(i,3), stress(i,4));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Float64" Name="Strain" NumberOfComponents="4" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e   %8.6e   %8.6e   %8.6e', strain(i,1), strain(i,2), strain(i,3), strain(i,4));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Float64" Name="Damage" NumberOfComponents="1" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e', damage(i,1));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Float64" Name="NLEquivalentStrain" NumberOfComponents="1" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e', NLEquivStrain(i,1));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Float64" Name="LS_psi" NumberOfComponents="1" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e', psi(i));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Float64" Name="LS_phi" NumberOfComponents="1" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, '  %8.6e', phi(i));
end 
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '</PointData>\n');
% fprintf(fid, '<CellData Scalars="" Vectors="" Tensors="" >\n');
% fprintf(fid, '</CellData>\n');
fprintf(fid, '</Piece>\n');


for i = 1:NELE
%     Trinodes = zeros(12,2);
    Tridisp = zeros(12,2);
    Tristress = zeros(12,4);
    Tridamage = zeros(12,1);
    Tristrain = zeros(12,4);
    TriNLEquivStrain = zeros(12,1);
    Tripsi = zeros(12,1);
    Triphi = zeros(12,1);
%     TriPorP = zeros(12,1);
%     TriVelocity = zeros(12,2);
    iElem = EnrichElems(i,1);
    N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
    N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
    N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
    N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
    NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element                 
    localD  = [N1*2-1 N1*2 N2*2-1 N2*2 N3*2-1 N3*2 N4*2-1 N4*2];    % Traditional index locations for displacement
%     localP  = [N1 N2 N3 N4];                                        % Traditional index locations for pressure
                                                                 % Next index location

    X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
    Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
    xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix
    if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
    end
    [Trinodes,Trinatural] = subDomain(PN,xyz);                        % Full Heaviside enrichment
    for j = 1:12
        xi = Trinatural(j,1); eta = Trinatural(j,2);                                % Gauss points
        N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...                     % Shape functions
                  (1+xi)*(1+eta);(1-xi)*(1+eta)];
        Tristress(j,1)=N'*stress([N1 N2 N3 N4]',1);
        Tristress(j,2)=N'*stress([N1 N2 N3 N4]',2);
        Tristress(j,3)=N'*stress([N1 N2 N3 N4]',3);
        Tristress(j,4)=N'*stress([N1 N2 N3 N4]',4);
        Tridamage(j,1)=N'*damage([N1 N2 N3 N4]',1);
        Tristrain(j,1)=N'*strain([N1 N2 N3 N4]',1);
        Tristrain(j,2)=N'*strain([N1 N2 N3 N4]',2);
        Tristrain(j,4)=N'*strain([N1 N2 N3 N4]',4);
        TriNLEquivStrain(j,1)=N'*NLEquivStrain([N1 N2 N3 N4]',1);
              
        Nenr = [];
        Nb = [];
        Nmatrix = [N(1,1)   0        N(2,1)   0        N(3,1)   0        N(4,1)   0;...
                   0        N(1,1)   0        N(2,1)   0        N(3,1)   0        N(4,1);];
            
        index = 1;
        iLocD   = 9;
        iLocP   = 5;
        
        psi_  = N(1)*PSI(N1)+N(2)*PSI(N2)+N(3)*PSI(N3)+N(4)*PSI(N4);        % Psi level set value at current point
        Tripsi(j,1) = psi_ ;
        Triphi(j,1) = N(1)*PHI(N1)+N(2)*PHI(N2)+N(3)*PHI(N3)+N(4)*PHI(N4);
        if ( abs(psi_) < 1e-6 )
            k = ceil(j/3);
            ind = (3*k-2):3*k;
            tricenter = sum(Trinatural(ind,:))./3;
            xi = tricenter(1,1); eta = tricenter(1,2);                      % triangle center
            Ncenter  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...                     % Shape functions
                      (1+xi)*(1+eta);(1-xi)*(1+eta)];
            psi_  = Ncenter(1)*PSI(N1)+Ncenter(2)*PSI(N2)+Ncenter(3)*PSI(N3)+Ncenter(4)*PSI(N4); 
        end
        
        for iN = 1:4
            if NN(iN,2) ~= 0
                Hgp = sign(psi_);                                                % Heaviside value at current gauss point
                Hi  = NN(iN,3);                                                 % Nodal Heaviside value
                H   = (Hgp-Hi)/2;                                               % Shifted Heaviside value
                Nenr(:,(2*index-1):(2*index)) = [ N(iN,1)*H    0;
                                                  0            N(iN,1)*H ];
%                 Hgp = abs(psi_);                                        % Heaviside value at current gauss point
%                 Hi  = abs(PSI(NN(iN,1)));                                         % Nodal Heaviside value
%                 D   = (Hgp-Hi);                                           % Shifted distance value
%                 Nb(:,index) = N(iN,1)*D;
                index = index+1;
                if ( j == 1 )
                    localD(iLocD:(iLocD+1)) = [2*NN(iN,2)-1 2*NN(iN,2)-0];
                    iLocD = iLocD+2;
%                     localP(iLocP) = NN(iN,2);
%                     iLocP = iLocP+1;
                end
            end
        end
        Nmat = [Nmatrix Nenr];    
        DSPD=DISPTD_(localD);
        Tridisp(j,:) = (Nmat*DSPD)';
%         ELEPORP=PRESTD_(localP);
%         TriPorP(j,:) = ([N' Nb]*ELEPORP)';
    end
	fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', 12,4);
    fprintf(fid, '<Points>\n');
    fprintf(fid, '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n');
    for i = 1:12
        fprintf(fid, ' %8.6e %8.6e %8.6e', Trinodes(i,1), Trinodes(i,2), 0);
    end
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '</Points>\n');
    fprintf(fid, '<Cells>\n');
    fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    for i = 1:4
        fprintf(fid, ' %d %d %d %d', i*3-3, i*3-2, i*3-1);
    end
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
    for i = 1:4
        fprintf(fid, ' %d', i*3);
    end
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
    for i = 1:4
        fprintf(fid, ' %d', 5);
    end
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, ' </Cells>\n');
    fprintf(fid, '<PointData Scalars="LS_psi LS_phi" Vectors="Displacement Stress Strain Damage NLEquivalentStrain">\n');
    fprintf(fid, '<DataArray type="Float64" Name="Displacement" NumberOfComponents="3" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e   %8.6e   %8.6e', Tridisp(i,1),Tridisp(i,2),0);
    end 
    fprintf(fid, ' </DataArray>\n');
%     fprintf(fid, '<DataArray type="Float64" Name="PorePressure" NumberOfComponents="1" format="ascii">\n');
%     for i = 1:12
%         fprintf(fid, '  %8.6e', TriPorP(i));
%     end 
%     fprintf(fid, ' </DataArray>\n');
%     fprintf(fid, '<DataArray type="Float64" Name="Velocity" NumberOfComponents="3" format="ascii">\n');
%     for i = 1:12
%         fprintf(fid, '  %8.6e   %8.6e   %8.6e', TriVelocity(i,1), TriVelocity(i,2), 0);
%     end 
%     fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="Stress" NumberOfComponents="4" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e   %8.6e   %8.6e   %8.6e', Tristress(i,1), Tristress(i,2), Tristress(i,3), Tristress(i,4));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="Strain" NumberOfComponents="4" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e   %8.6e   %8.6e   %8.6e', Tristrain(i,1), Tristrain(i,2), Tristrain(i,3), Tristrain(i,4));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="Damage" NumberOfComponents="1" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e', Tridamage(i,1));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="NLEquivalentStrain" NumberOfComponents="1" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e', TriNLEquivStrain(i,1));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="LS_psi" NumberOfComponents="1" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e', Tripsi(i));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '<DataArray type="Float64" Name="LS_phi" NumberOfComponents="1" format="ascii">\n');
    for i = 1:12
        fprintf(fid, '  %8.6e', Triphi(i));
    end 
    fprintf(fid, ' </DataArray>\n');
    fprintf(fid, '</PointData>\n');
    fprintf(fid, '</Piece>\n');
end
fprintf(fid, '</UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');
fclose(fid);

% I_frac = find( connec_frac(:,2) == 1 );
% CMOD = statev_frac{I_frac}{1}.width;
% ind=size(xyz_frac,1);
% cracklength = xyz_frac(ind,4);
% filename = strcat('fracPressure.txt');
% fid = fopen(filename, 'a');
% fprintf(fid, '  %8.6e   %8.6e   %8.6e', HydroPressure(1),CMOD,cracklength); 
% fprintf(fid, '\n');
% fclose(fid);

ind=find(EXTDISP(:,2)~=0);
Displacement_node=EXTDISP(ind,1);
FORCE_ = full(FORCE);
filename = strcat(file_name,'Reaction.txt');
fid = fopen(filename, 'a');
for i=1:size(Displacement_node,1)
    fprintf(fid, '   %8.6e   %8.6e', DISPTD_(Displacement_node(i)),FORCE_(Displacement_node(i))); 
end
fprintf(fid, '\n');
fclose(fid);

end


function [Trinodes,Trinatural] = subDomain(psi,xyz)
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


Trinodes = zeros(12,2);
Trinatural = zeros(12,2);

for e = 1:size(tri,1)
    coord = node(tri(e,:),:);
    xyzl  = xyz(tri(e,:),:);
    index = 3*e-2:e*3;
    Trinodes(index,:) = xyzl ;
    Trinatural(index,:) = coord ;
end

end

