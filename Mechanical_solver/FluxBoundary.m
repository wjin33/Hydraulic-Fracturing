% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function [EXTFlUX]=FluxBoundary(Bsurface,surf,applied_flux,injection_rate)
% This function creates the global force vector from the user supplied
% input file.

global CONNEC XYZ NODES

n = size(surf,2);      % Number of boundary surface with applied stress

m = size(Bsurface,1);  % Number of surface defined in the input file

% Create the global force matrix from distributed flux, we do not use to for the hydraulic pressurization case
globalF = zeros(3*size(XYZ,1),1);                                           
for iForce = 1:m
    for iSurf = 1:n
        if (strcmp(Bsurface{iForce,1},surf(iSurf))>0)
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

                force = applied_flux(iSurf)*edgelength;
            
                globalF(edge_coordinates(1,1),1) = globalF(edge_coordinates(1,1),1) + force/2;   % Share of node from Elemental flux flow into the domain                    des
                globalF(edge_coordinates(2,1),1) = globalF(edge_coordinates(2,1),1) + force/2; 
            end
            break;
        end
    end
end

index = find(globalF(:,1)~=0);
EXTFlUX = zeros(size(index,1), 2);

EXTFlUX(:,1) = index;
EXTFlUX(:,2) = globalF(index,1);

% Apply the fluid injection rate over time, need to be implemented other if
% injection rate changes over time

EXTFlUX(1,:) = [max(max(NODES))+1 injection_rate];

  
end