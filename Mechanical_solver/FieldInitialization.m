% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function [ NEQ ] = FieldInitialization()
%%
global PREFTD NODES xyz_frac

% NEQ=[0 0];
NEQ = 2*max(max(NODES));             % Number of degrees of freedom for displacement
% NEQ(2) = size(xyz_frac,1);                               % Number of degrees of freedom for pressure
PREFTD  = sparse(NEQ(1),1);                              % The displacement as well as the pore pressure is initially zero;
% HydroPressure  = ones(NEQ(2),1);
% %% 

heavDOF = nnz(NODES(:,2));                                                % Define the number of Heavi DOF
if heavDOF > 0, heaviNodes;     end

end

function heaviNodes
% This function assigns nodes enriched with the Heaviside function as
% either above (+1) or below (-1) the crack.

global NODES PSI

for iNode = 1:size(NODES,1)
    if NODES(iNode,2) ~= 0
        NODES(iNode,3) = sign(PSI(iNode));
    end
end
end