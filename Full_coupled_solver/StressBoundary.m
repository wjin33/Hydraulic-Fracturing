% Written By: Wencheng Jin, Idaho National Laboratory (2020)
% Website: https://sites.google.com/view/wenchengjin/software
% Email: wencheng.jin@inl.gov

function [EXTFORCE]=StressBoundary(Bsurface,surf,stressT)
% This function creates the global force vector from the user supplied input file.

global CONNEC XYZ

n = size(surf,2);      % Number of boundary surface with applied stress

m = size(Bsurface,1);  % Number of surface defined in the input file

globalF = zeros(3*size(XYZ,1),1);                                              % Create the global force vector

for iForce = 1:m
    for iSurf = 1:n
        if (strcmp(Bsurface{iForce,1},surf(iSurf))>0)
            surface = Bsurface{iForce,2};
            num_edge = size(surface ,2)/2;
            for iedge = 1:num_edge
                Nodes  = CONNEC(surface(iedge*2-1),2:5);                        % Nodes for current element
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

                force = stressT(iSurf,:)*edgelength;
        
                globalF(3*edge_coordinates(1,1)-2,1) = globalF(3*edge_coordinates(1,1)-2,1) + force(1)/2;                      % Elemental force in x-direction for edge nodes
                globalF(3*edge_coordinates(1,1)-1,1) = globalF(3*edge_coordinates(1,1)-1,1) + force(2)/2;                      % Elemental force in y-direction for edge nodes
                globalF(3*edge_coordinates(2,1)-2,1) = globalF(3*edge_coordinates(2,1)-2,1) + force(1)/2;                      % Elemental force in x-direction for edge nodes
                globalF(3*edge_coordinates(2,1)-1,1) = globalF(3*edge_coordinates(2,1)-1,1) + force(2)/2;   
            end
%             break;
        end
    end
end

index = find(globalF(:,1)~=0);
EXTFORCE = zeros(size(index,1), 2);

EXTFORCE(:,1) = index;
EXTFORCE(:,2) = globalF(index,1);
  
end