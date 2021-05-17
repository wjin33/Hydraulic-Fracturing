% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function [EXTDISP,EXTFlUX,EXTPRESSURE]=FluxBoundary(Bsurface,surf,applied_flux,injection_rate,EXTDISP,EXTPRESSURE)
% This function creates the global force vector from the user supplied
% input file.

global CONNEC XYZ NODES

n = size(surf,2);      % Number of boundary surface with applied stress

m = size(Bsurface,1);  % Number of surface defined in the input file

% Create the global force matrix from distributed flux, we do not use to for the hydraulic pressurization case
globalF = zeros(3*max(max(NODES)),1);                                           
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
%                 edgelength = sqrt((edge_coordinates(2,2)-edge_coordinates(1,2))^2 + (edge_coordinates(2,3)-edge_coordinates(1,3))^2);
                ind = [NODES(edge_coordinates(1,1),2); NODES(edge_coordinates(2,1),2)];            
                globalF(3*edge_coordinates(1,1),1) = injection_rate/2;   % Share of node from Elemental flux flow into the domain                    des
                globalF(3*edge_coordinates(2,1),1) = injection_rate/2; 
            end
            break;
        end
    end
end

index = find(globalF(:,1)~=0);
EXTFlUX = zeros(size(index,1), 2);
EXTFlUX(:,1) = index;
EXTFlUX(:,2) = globalF(index,1);

EXTDISP = [EXTDISP; ind(1)*3-2, 0; ind(2)*3-2 0;];

EXTPRESSURE = [EXTPRESSURE;ind(1)*3, 0; ind(2)*3 0; ];

return;

%Apply node flux, fluid injection for hydraulic fracturing

% n = size(set,2);      % Number of node sets with applied pressure, drained condition                                     
% m = size(BNoset,1);     % Number of node sets defined in input file
% Tempflux = [];
% for iPressure = 1:m
%     for iSet = 1:n
%         if (strcmp(BNoset{iPressure,1},set(iSet))>0)
%             nodes = BNoset{iPressure,2};
%             temp = ones(size(nodes));
%             Tempflux = [Tempflux; [3*nodes; temp]'];
%         end
%     end
% end
% 
% m = sum(Tempflux(:,2));
% Tempflux(:,2) = ones(m,1).*injection_rate/m;

% index = find(NODES(:,2)~=0);
% m = size(index,1);
% Tempflux=zeros(m,2);
% Tempflux(:,1) = index.*3;
% Tempflux(:,2) =  ones(m,1).*(injection_rate/m);
% EXTFlUX =[EXTFlUX;Tempflux];

for iForce = 1:m
    for iSurf = 1:n
        if (strcmp(Bsurface{iForce,1},surf(iSurf))>0)
            surface = Bsurface{iForce,2};
            num_edge = size(surface ,2)/2;
            for iedge = 1:num_edge
                Nodes  = CONNEC(surface(iedge*2-1),2:5);                                                  % Nodes for current element
                ind = [NODES(Nodes',2)];
                globalF(3*Nodes(1,1),1) =  injection_rate/4;   % Share of node from Elemental flux flow into the domain                    des
                globalF(3*Nodes(1,2),1) =  injection_rate/4;
                globalF(3*Nodes(1,3),1) =  injection_rate/4;
                globalF(3*Nodes(1,4),1) =  injection_rate/4;
            end
            break;
        end
    end
end


% TempDisp = zeros(m,2);
% TempDisp(:,1) = NODES(Tempflux(:,1)./3,2)*3-2;
% EXTDISP = [EXTDISP;TempDisp];
EXTDISP =EXTDISP;
end