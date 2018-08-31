% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function InitialBoundary(Bsurface)
% This function creates the global force vector from the user supplied
% input file.

global Confinement CONNEC XYZ FORCE

% globalF = sparse(NEQ,1);                                                    % Create the global force vector
m = size(Confinement,2);

% Create the global force matrix from distributed loading

for iForce = 1:m
    surface = Bsurface{iForce,2};
    num_edge = size(surface ,2)/2;
    for iedge = 1:num_edge
        Nodes  = CONNEC(surface(iedge*2-1),2:5);                                                  % Nodes for current element
        switch surface(iedge*2)
            case 1
                edge_coordinates  = XYZ([Nodes(1) Nodes(2)]',:);
            case 2
                edge_coordinates  = XYZ([Nodes(2) Nodes(3)]',:);
            case 3
                edge_coordinates  = XYZ([Nodes(3) Nodes(4)]',:);
            case 4
                edge_coordinates  = XYZ([Nodes(4) Nodes(1)]',:);
        end
        edgelength = sqrt((edge_coordinates(2,2)-edge_coordinates(1,2))^2 + (edge_coordinates(2,3)-edge_coordinates(1,3))^2);
        vector = edge_coordinates(2,2:3) - edge_coordinates(1,2:3);
        if sum(vector ~= 0) == 2
            v_perp(1) = 1;
            v_perp(2) = -vector(1)/(vector(2));
            v_perp = v_perp./sqrt(v_perp(2)^2+ v_perp(1)^2);
        else
            if vector(1) ~= 0
               v_perp(1) = 0;
               v_perp(2) = 1;
            else
               v_perp(1) = 1;
               v_perp(2) = 0;
            end
        end
        
        if (v_perp(2)*vector(1)-v_perp(1)*vector(2) < 0 )
            v_perp = -v_perp;
        end
        force = Confinement(iForce)*edgelength*v_perp;
        
        FORCE(2*edge_coordinates(1,1)-1,1) = FORCE(2*edge_coordinates(1,1)-1,1) + force(1)/2;                      % Elemental force in x-direction for edge nodes
        FORCE(2*edge_coordinates(1,1),1)   = FORCE(2*edge_coordinates(1,1),1) + force(2)/2;                      % Elemental force in y-direction for edge nodes
        FORCE(2*edge_coordinates(2,1)-1,1) = FORCE(2*edge_coordinates(2,1)-1,1) + force(1)/2;                      % Elemental force in x-direction for edge nodes
        FORCE(2*edge_coordinates(2,1),1)   = FORCE(2*edge_coordinates(2,1),1) + force(2)/2;                      % Elemental force in y-direction for edge nodes 
   
    end
end
   
end