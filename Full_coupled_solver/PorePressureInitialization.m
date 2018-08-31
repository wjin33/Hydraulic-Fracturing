function [ NEQ ] = PorePressureInitialization(EXTPRESSURE)
%%
global PREPTD NODES xyz_frac
%Initial condition of pore pressur distribution

NEQ = max(xyz_frac(:,1));
if NEQ == 0
    NEQ = max(max(NODES));                                                % Total number of degrees of freedom
end
PREPTD = zeros(NEQ,1);
NPRE=size(EXTPRESSURE,1);
if NPRE~=0
    FIXEDDOF=EXTPRESSURE(:,1);
    PREPTD(FIXEDDOF) = EXTPRESSURE(:,2);
end

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