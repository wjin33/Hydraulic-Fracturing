% Written By: Wencheng Jin, Idaho National Laboratory (2020)
% Website: https://sites.google.com/view/wenchengjin/software
% Email: wencheng.jin@inl.gov

function [enrElem] = enrElem
% This function finds the enriched elements.

global CONNEC NODES

    % Find the new elements containing crack body domain
    I = find(NODES(:,2) > 0);                                         % Nodes containing crack tip enrichment
    [~,~,c1]    = intersect(I,CONNEC(:,2)');                                % Find elements with crack tip enrichment
    [~,~,c2]    = intersect(I,CONNEC(:,3)');                                % Find elements with crack tip enrichment
    [~,~,c3]    = intersect(I,CONNEC(:,4)');                                % Find elements with crack tip enrichment
    [~,~,c4]    = intersect(I,CONNEC(:,5)');                                % Find elements with crack tip enrichment
    enrElem    = [];                                                       % Elements fully enriched
    tipelement = [];
    candidate = unique([c1; c2; c3; c4;]);
    for i = 1:size(candidate, 1)
        iElem = candidate(i);
        N1  = CONNEC(iElem,2);                                              % Node 1 for current element
        N2  = CONNEC(iElem,3);                                              % Node 2 for current element
        N3  = CONNEC(iElem,4);                                              % Node 3 for current element
        N4  = CONNEC(iElem,5);                                              % Node 4 for current element
        nodes = [N1 N2 N3 N4]';
        NN  = NODES(nodes,:);                                               % Nodal data for current element
        SIGN = abs(sum(NN(:,3)));
        HEN = nnz(NN(:,2));
        if ( HEN == 4 && SIGN ~=4 )
            enrElem = [enrElem;iElem];
        end
%         if ( HEN == 2 )
%             if sum(NN(:,3)) == 0
%                 tipelement = iElem;
% %                 break;
%             end
%         end
    end
    enrElem = [enrElem; tipelement];
end