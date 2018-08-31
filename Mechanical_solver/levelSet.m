% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu
function []=levelSet(update)
% This function creates the ZETA level set function representing the
% materials within the domain. This function creates the CHI level set
% function representing voids within the domain. This function creates the
% PHI and PSI level set functions needed to represent the discontinuity
% created by the crack. The elements with Heaviside and crack tip
% enrichments are determined and node numbers are assigned to the nodes
% which have enriched degrees of freedom.

global   CRACK CONNEC NODES XYZ PHI PSI 

nElem  = size(CONNEC,1);                                            % Number of elements in domain
nNode  = size(NODES,1);                                             % Number of nodes in domain
nPt    = size(CRACK,1);                                             % Number of data points defining crack

Xmin=min(XYZ(:,2));
Xmax=max(XYZ(:,2));
Ymin=min(XYZ(:,3));
Ymax=max(XYZ(:,3));

% Create the level set functions (PHI and PSI) defining the crack
if isempty(CRACK) == 0
    nCT = 2;                                                                % Default number of crack tips
    if     (CRACK(1,1)   <= Xmin) || (CRACK(1,1)   >= Xmax)                 % Check for edge crack
        nCT = nCT-1;
    elseif (CRACK(nPt,1) <= Xmin) || (CRACK(nPt,1) >= Xmax)
        nCT = nCT-1;
    elseif (CRACK(1,2)   <= Ymin) || (CRACK(1,2)   >= Ymax)
        nCT = nCT-1;
    elseif (CRACK(nPt,2) <= Ymin) || (CRACK(nPt,2) >= Ymax)
        nCT = nCT-1;
    end

    PHI = sparse(nNode,nCT);                                                % Initialize phi
    if (~update==1)                                                            % Initialize psi
        PSI = sparse(nNode,1);
    end 
    
    segcenter=sum(CRACK)./nPt;
    crack_len=0;
    % Find the center of crack segment
    for i=2:nPt
        crack_len = crack_len + sqrt((CRACK(i-1,1)-CRACK(i,1))^2+(CRACK(i-1,2)-CRACK(i,2))^2);
    end
                                                    % search radius
    Distance2Tip = sqrt( (XYZ(:,2)-CRACK(nPt,1)).^2 + (XYZ(:,3)-CRACK(nPt,2)).^2);
    [~,I] = min(Distance2Tip);
    [~,~,c1] = intersect(I,CONNEC(:,2)');
    TipElementNodes = CONNEC(c1(1),:);
    elementSize = sqrt( (XYZ(TipElementNodes(2),2)-XYZ(TipElementNodes(4),2))^2 + (XYZ(TipElementNodes(2),3)-XYZ(TipElementNodes(4),3)).^2);
    radius = crack_len+3*elementSize;
    
%     radius =  crack_len;
       
    dist = zeros(1,nNode);                                                  % Initialize distance vector
    xCT    = segcenter(1,1);                                                % X-coordinate of crack center
    yCT    = segcenter(1,2);                                                % Y-coordinate of crack cente  
    for iN = 1:nNode
        Xn       = XYZ(iN,2);                                               % X-coordinate for the current node
        Yn       = XYZ(iN,3);                                               % Y-coordinate for the current node
        dist(iN) = sqrt((Xn-xCT)^2+(Yn-yCT)^2);                             % Store radius value
    end

    temp   = dist-radius;                                                   % Determine whether or not the node is outside the search radius
    domain = find(temp <= 0);                                               % Find nodes within search radius
        
    % Compute the PHI level set functions for the main crack tip(s)
    for iNode = 1:length(domain)
        cNode = domain(iNode);                                              % Current node within search radius
        x     = XYZ(cNode,2);                                               % X-coordinate for the current node
        y     = XYZ(cNode,3);                                               % Y-coordinate for the current node
        
        if nCT == 1
            % Define phi for first crack tip
            disc = CRACK(nPt,:)-CRACK(nPt-1,:);                             % Horizontal and vertical distances for current crack segment
            t    = 1/norm(disc)*disc;                                       % Tangent to current crack segment
            phi  = ([x y]-CRACK(nPt,:))*t';
            if phi == 0, phi = 1e-6; end
            % Define phi and psi at nodes within narrow band
            dist = zeros(nPt-1,1); sine = dist;
            for iSeg = 1:(nPt-1)
                x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        
                m = (y2-y1)/(x2-x1);
                if isinf(m) == 1
                    dist(iSeg) = (abs(x-x1)+abs(x-x2))/2;
                    vec = [x y]-[x1 y1];
                    sine(iSeg) = sign(t(1)*vec(2)-t(2)*vec(1));
                    continue;
                end
                b = y1-m*x1;
                        
                xo = (m*y+x-m*b)/(m^2+1);
                yo = (m^2*y+m*x+b)/(m^2+1);
                        
                if iSeg ~= 1
                    if ( xo - x1) < -1e-6
                        xo = x1;
                        yo = y1;
                    end
                end
                if iSeg ~= (nPt-1)
                    if ( xo-x2 > 1e-6) 
                        xo = x2;
                        yo = y2;
                    end
                end      
                dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                sine(iSeg) = sign(x2-x1)*sign(y-yo);
                if sign(x2-x1) == 0
                    sine(iSeg) =sign(y-yo);
                end
            end
            dMin = min(abs(dist));
            ind  = find(abs(dist) == dMin);
            if length(ind) >= 2, ind = ind(1); end
            psi = dist(ind)*sine(ind);
            if psi == 0, psi = 1e-6; end    
            PHI(cNode,1) = phi(1);
            PSI(cNode,1) = psi;
            
            
        elseif ( nCT == 2 )
            % Define phi for first crack tip
            disc = CRACK(nPt,:)-CRACK(nPt-1,:);                             % Horizontal and vertical distances for current crack segment
            t    = 1/norm(disc)*disc;                                       % Tangent to current crack segment
            phi  = ([x y]-CRACK(nPt,:))*t';
            if phi == 0, phi = 1e-6; end
            % Define phi and psi at nodes within narrow band
            dist = zeros(nPt-1,1); sine = dist;
            for iSeg = 1:(nPt-1)
                x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        
                m = (y2-y1)/(x2-x1);
                if isinf(m) == 1, m = 1e6; end
                b = y1-m*x1;
                        
                xo = (m*y+x-m*b)/(m^2+1);
                yo = (m^2*y+m*x+b)/(m^2+1);
                        
                if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                sine(iSeg) = sign(x2-x1)*sign(y-yo);
            end
            dMin = min(abs(dist));
            ind  = find(abs(dist) == dMin);
            if length(ind) == 2, ind(2) = []; end
            psi = dist(ind)*sine(ind);
            if psi == 0, psi = 1e-6; end    
            PHI(cNode,1) = phi(1);
            PSI(cNode,1) = psi;
            
            disc = CRACK(1,:)-CRACK(2,:);                               % Horizontal and vertical distances for current crack segment
            t    = 1/norm(disc)*disc;                                   % Tangent to current crack segment
            phi(2)  = ([x y]-CRACK(1,:))*t';
            if phi(2) == 0, phi(2) = 1e-6; end
            if PHI(cNode,1) < phi(2)
                PHI(cNode,1) = phi(2);
            end
        end    
    end
end
    

    
    % Define Heaviside enriched nodes
    hNodes   = [];
    I        = find(PSI ~= 0);                                              % Nodes with defined psi
    [~,~,c1] = intersect(I,CONNEC(:,2)');                                   % Find elements with defined psi
    [~,~,c2] = intersect(I,CONNEC(:,3)');                                   % Find elements with defined psi
    [~,~,c3] = intersect(I,CONNEC(:,4)');                                   % Find elements with defined psi
    [~,~,c4] = intersect(I,CONNEC(:,5)');                                   % Find elements with defined psi
    hElem    = unique([c1' c2' c3' c4']);                                       % Candidate elements for Heaviside enrichment
    
    for iElem = 1:length(hElem)
        cElem = hElem(iElem);
        psiE  = PSI(CONNEC(cElem,2:5));
        if nnz(psiE) == 4
            if (max(psiE) == 1e-6) || (min(psiE) == 1e-6)
                for iN = 1:4
                    gN = CONNEC(cElem,iN+1);
                    if psiE(iN) == 1e-6
                        if PHI(gN) <= 0
                            if NODES(gN,2) == 0
                                hNodes = [hNodes gN];
                            end
                        end
                    end
                end
            elseif max(psiE)*min(psiE) < 0
                for iN = 1:4
                    gN = CONNEC(cElem,iN+1);
                    if PHI(gN) <= 0
                        if NODES(gN,2) == 0
                            hNodes = [hNodes gN];
                        end
                    end
                end
            end
        end
    end
    
    % determine the previous enriched maximum node number 
    if update
        incNode = 1+max(NODES(:,2))-nNode;
    else
        incNode = 1;
    end
    
    % Preliminaries to numbering enriched nodes
    nNode   = nNode+incNode;
    hNodes  = unique(hNodes);
    
    % Number the Heaviside nodes
    for i = 1:length(hNodes)
        NODES(hNodes(i),2) = nNode;
        nNode = nNode+1;
    end

end