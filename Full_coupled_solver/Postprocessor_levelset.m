% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function Postprocessor_levelset()
% This function plots the PHI, PSI and ZETA level set functions.

global  PHI PSI XYZ CONNEC 
%Recover the state varibals at the nodes from Gauss Points

psi = full(PSI);
phi = full(PHI);
npoints=size(XYZ,1);
ncells=size(CONNEC,1);


filename = 'doubleTip.vtu';
fid = fopen(filename, 'w');
% fprintf(fid, '<!-- TimeStep %8.6e -->\n',TIME);
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '<UnstructuredGrid>\n');
fprintf(fid, '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n', npoints,ncells);
fprintf(fid, '<Points>\n');
fprintf(fid, '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n');
for i = 1:npoints
    fprintf(fid, ' %8.6e %8.6e %8.6e', XYZ(i,2), XYZ(i,3), 0);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '</Points>\n');
fprintf(fid, '<Cells>\n');
fprintf(fid, '<DataArray type="Int32" Name="connectivity" format="ascii">\n');
for i = 1:ncells
    iEle = i;
    fprintf(fid, ' %d %d %d %d', CONNEC(iEle,2)-1, CONNEC(iEle,3)-1, CONNEC(iEle,4)-1,CONNEC(iEle,5)-1);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="Int32" Name="offsets" format="ascii">\n');
for i = 1:(ncells)
    fprintf(fid, ' %d', i*4);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, '<DataArray type="UInt8" Name="types" format="ascii">\n');
for i = 1:(ncells)
    fprintf(fid, ' %d', 9);
end
fprintf(fid, ' </DataArray>\n');
fprintf(fid, ' </Cells>\n');
% fprintf(fid, '<PointData Scalars="LS_psi LS_phi " Vectors="DisplacementVector " Tensors="StressTensor DamageTensor " >\n');
% fprintf(fid, '<PointData Scalars="LS_psi LS_phi" Vectors="DisplacementVector StressTensor DamageTensor EquivalentStrain">\n');
fprintf(fid, '<PointData Scalars="LS_psi LS_phi" Vectors="Displacement Stress Strain Damage NLEquivalentStrain">\n');
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
fprintf(fid, '</UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');
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

