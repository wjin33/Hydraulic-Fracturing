% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function [EXTPRESSURE]=PressureBoundary(BNoset,set,pressureB)

global NODES

n = size(set,2);      % Number of node sets with applied pressure, drained condition
                                         
m = size(BNoset,1);     % Number of node sets defined in input file

EXTPRESSURE = [];
for iPressure = 1:m
    for iSet = 1:n
        if (strcmp(BNoset{iPressure,1},set(iSet))>0)
            nodes = BNoset{iPressure,2};
            temp = pressureB(iSet)*ones(size(nodes));
            EXTPRESSURE = [EXTPRESSURE; [3*nodes; temp]'];
        end
    end
end


% EXTPRESSURE(1,:) = [3*max(max(NODES))+1 10];
  
end