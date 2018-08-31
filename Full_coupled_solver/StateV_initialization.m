% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu


function StateV_initialization(PROP)	% Initialize elastic history dependent varibles
%%
  global STATEV CONNEC NODES ELECENTER XYZ
  %
  eqeps_1t = PROP.eqeps_1t;
  eqeps_2t = PROP.eqeps_2t;
  NE = size(CONNEC, 1);
  STATEV=cell(1,NE);
  ELECENTER = zeros(NE,2);
  % history.numbEle = [0];
  history.natural_coodinates = [0; 0];
  history.gauss_weight = 0;
  history.true_coodinates = [0; 0];
  history.sigma = [0; 0; 0; 0;];
  history.strain = [0; 0; 0; 0;];
  history.damage = [0; 0];
  history.kappa = [eqeps_1t; eqeps_2t];
  history.EquivStrain = [0; 0];
  history.NLEquivStrain = [0; 0];
  history.nonlocalTable = [];
  history.volume = 0;
  history.fluidVelocity = [0; 0];
  %

  %
  for iElem = 1:NE
      N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
      N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
      N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
      N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
      nodes = [N1 N2 N3 N4]';
      NN  = NODES(nodes,:);                                          % Nodal data for current element
      ELECENTER(iElem,:) = sum(XYZ(nodes,2:3))./4;  
      HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
      
      if HEN == 4  
          STATEV{iElem}={history, history, history, history,history, history, history, history,history, history, history, history};
%       elseif HEN == 2
%           STATEV{iElem}={history, history, history, history, history, history, history, history, history};
      else
          STATEV{iElem}={history, history, history, history};
      end
  end
  
  %initialize natrual and true coordinates
  GPInitialization;
  
end

function  GPInitialization()
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

  global CONNEC NODES PSI XYZ STATEV connec_frac xyz_frac statev_frac CRACK PROP

  fracxy = [];
  enrElems = [];

  for iElem = 1:size(CONNEC,1)
      N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
      N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
      N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
      N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
      NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
      HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    
      if (HEN == 0)                                                           % Unenriched nodes
          % Traditional element
          X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2); % Nodal x-coordinates
          Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3); % Nodal y-coordinates     
          xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix   
        
          [gp,gw] = gauss(2,'QUAD');
          for ig = 1:length(gp)
              xi = gp(ig,1); eta = gp(ig,2);                                % Gauss points
              STATEV{iElem}{ig}.natural_coodinates = [xi;eta];
              N = 1/4*[(1-xi)*(1-eta)   (1+xi)*(1-eta)  (1+xi)*(1+eta)  (1-xi)*(1+eta)];          % Shape functions
              STATEV{iElem}{ig}.true_coodinates = [N*xyz]';
              STATEV{iElem}{ig}.gauss_weight = gw(ig,1);
          end   
        
      elseif HEN > 0                                                                 % Enriched element
        
          X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
          Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
          xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix

          if HEN == 4                                                         % Fully enriched element
              if numel(PSI) == 0, PN = [0 0 0 0]; else
                  PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
              end
              [gp,gw,~,fracnode] = subDomain(3,PN,xyz);                         % Full Heaviside enrichment
              fracxy = [fracxy; fracnode];
              enrElems = [enrElems iElem];
          else                                                                % Partially enriched element
              [gp,gw] = gauss(2,'QUAD');
          end

          for ig = 1:length(gp)
              xi = gp(ig,1); eta = gp(ig,2);                                    % Gauss points\
            
              STATEV{iElem}{ig}.natural_coodinates = [xi;eta];
              N = 1/4*[(1-xi)*(1-eta)   (1+xi)*(1-eta)  (1+xi)*(1+eta)  (1-xi)*(1+eta)];          % Shape functions
              STATEV{iElem}{ig}.true_coodinates = [N*xyz]';
              STATEV{iElem}{ig}.gauss_weight = gw(ig,1);
          end
      end
  end
  
%   return;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  frac.natural_coodinates = 0;
  frac.gauss_weight = 0;
  frac.true_coodinates = [0; 0];
  frac.natural_coordinate_ele = [0;0];
  frac.width = 0.0;         %mm, tentatively, should be updated from mechanical calculation
  frac.fluidVelocity = 0;
  frac.seglength = 0;
  frac.jump = [1e6;1e6];
  frac.GI = 0;                  %N/m
  frac.GII = 0;                  %N/m
  frac.sigmaMax = 0;              %N/m^2  Pa
  frac.tauMax = 0;              %N/m^2  Pa
  frac.lambdaN = 0;
  frac.lambdaT = 0;
  frac.alpha = 0;
  frac.beta = 0;
  frac.m = 0;
  frac.n = 0;
  frac.deltaN = 0;
  frac.deltaT = 0;
  frac.PenaltyStiffness = 0;
  frac.dGnt = 0;
  frac.dGtn = 0;
  frac.deltaN_conj = 0;
  frac.deltaT_conj = 0;
  frac.GammaN = 0;
  frac.GammaT = 0;
  segment_number = size(enrElems,2);
  statev_frac = cell(1,segment_number);
  connec_frac = zeros(segment_number,3);
  xyz_frac = zeros(segment_number+1,4);
  maxNodeNum = max(NODES(:,2));
  xyz_frac(:,1) = [(0*maxNodeNum+1) : (0*maxNodeNum + segment_number+1)]';
  UnorderedNodes = unique(fracxy,'rows');
  backup_UnorderedNodes =UnorderedNodes;
  orderedNodes = [];
  
  k = dsearchn(UnorderedNodes,delaunayn(UnorderedNodes),CRACK(1,:));
  orderedNodes = UnorderedNodes(k,:);
  UnorderedNodes = setdiff(UnorderedNodes,orderedNodes,'rows');
  while (size(UnorderedNodes,1)>2)
      lambda = size(orderedNodes,1);
      k = dsearchn(UnorderedNodes,delaunayn(UnorderedNodes),orderedNodes(lambda,:));
      orderedNodes =[ orderedNodes;  UnorderedNodes(k,:)];
      UnorderedNodes = setdiff(UnorderedNodes,UnorderedNodes(k,:),'rows');
  end
  lambda = size(orderedNodes,1);
  dist1 = sqrt( sum( (UnorderedNodes(1,:)-orderedNodes(lambda,:)).^2 ) );
  dist2 = sqrt( sum( (UnorderedNodes(2,:)-orderedNodes(lambda,:)).^2 ) );
  if dist1 < dist2
      orderedNodes = [orderedNodes; UnorderedNodes(1,:) ; UnorderedNodes(2,:) ];
  else
      orderedNodes = [orderedNodes; UnorderedNodes(2,:) ; UnorderedNodes(1,:) ];
  end
%       
%   nPt = size(CRACK,1);
%   for iSeg = 1:(nPt-1)
%       x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
%       x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
%       vectorA = [CRACK(iSeg+1,1)-CRACK(iSeg,1), CRACK(iSeg+1,2)-CRACK(iSeg,2)];
%       TempNodes = [];
%       lambda = [];
%       for inode = 1:size(UnorderedNodes,1)
%           vectorB = [UnorderedNodes(inode,1)-CRACK(iSeg,1), UnorderedNodes(inode,2)-CRACK(iSeg,2)];
%           if vectorA(1) == 0
%               if vectorB(1) == 0
%                   lambday = vectorB(2)/vectorA(2);
%                   TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                   lambda = [lambda; lambday];
%               end
%           elseif vectorA(2) == 0
%               if vectorB(2) == 0
%                   lambdax = vectorB(1)/vectorA(1);
%                   TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                   lambda = [lambda; lambdax];
%               end
%           else     
%             lambdax = vectorB(1)/vectorA(1);
%             lambday = vectorB(2)/vectorA(2); 
%             if abs( lambdax - lambday ) <= 1E-4
%                 TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                 lambda = [lambda; lambday];
%             end
%           end
%       end
%       if size(lambda,1)~=0
%           [~,I] = sort(lambda);
%           UnorderedNodes = setdiff(UnorderedNodes,TempNodes,'rows');
%           orderedNodes = [orderedNodes ; TempNodes(I,:)];
%       end
%   end
  if size(orderedNodes,1) == size(xyz_frac,1)
      xyz_frac(:,2:3)  = orderedNodes;
  else
      xyz_frac(:,2:3)  = backup_UnorderedNodes;
  end
%  xyz_frac(:,2:3)  = UnorderedNodes;
  for inode = 1 : segment_number
      xyz_frac(inode+1,4) = xyz_frac(inode,4) + sqrt(sum((xyz_frac(inode+1,2:3)-xyz_frac(inode,2:3)).^2));
  end
  
  for i = 1:segment_number
      xynode = fracxy(2*i-1:2*i, :);
      Ia = find(xyz_frac(:,2) == xynode(1,1));
      Ib = find(xyz_frac(:,3) == xynode(1,2));
      node1 = intersect(Ia,Ib);
      Ia = find(xyz_frac(:,2) == xynode(2,1));
      Ib = find(xyz_frac(:,3) == xynode(2,2));
      node2 = intersect(Ia,Ib);
      connec_frac(i,:) = [enrElems(i) sort([node1 node2])];
      slope = xynode(1,:)-xynode(2,:);
      theta = atan(slope(2)/slope(1));
      frac.GI = sqrt((PROP.G2*cos(theta))^2 + (PROP.G1*sin(theta))^2);                  %N/m
      frac.GII = frac.GI;                  %N/m
      frac.sigmaMax = sqrt((PROP.sigmaMax2*cos(theta))^2 + (PROP.sigmaMax1*sin(theta))^2);              %N/m^2  Pa
      frac.tauMax = frac.sigmaMax;              %N/m^2  Pa
      frac.lambdaN = 0.01;
      frac.lambdaT = 0.01;
      frac.alpha = 4;
      frac.beta = 4;
      frac.m = frac.alpha*(frac.alpha-1)*frac.lambdaN^2/(1-frac.alpha*frac.lambdaN^2);
      frac.n = frac.beta*(frac.beta-1)*frac.lambdaT^2/(1-frac.beta*frac.lambdaT^2);
      frac.deltaN = frac.GI/frac.sigmaMax*frac.alpha*frac.lambdaN*(1-frac.lambdaN)^(frac.alpha-1)*(frac.alpha/frac.m+1)*(frac.lambdaN*frac.alpha/frac.m+1)^(frac.m-1);
      frac.deltaT = frac.GII/frac.tauMax*frac.alpha*frac.lambdaT*(1-frac.lambdaT)^(frac.beta-1)*(frac.beta/frac.n+1)*(frac.lambdaT*frac.beta/frac.n+1)^(frac.n-1);
      frac.PenaltyStiffness = 1e8*frac.sigmaMax/frac.deltaN;
      frac.dGnt = 0;
      frac.dGtn = 0;
      frac.deltaN_conj = frac.deltaN-frac.deltaN*(frac.dGnt/frac.GI)^(1/frac.alpha);
      frac.deltaT_conj = frac.deltaT-frac.deltaT*(frac.dGtn/frac.GI)^(1/frac.beta);
      frac.GammaN = -frac.GI*(frac.alpha/frac.m)^frac.m;
      frac.GammaT = (frac.beta/frac.n)^frac.n ;
      frac.seglength = sqrt(sum((xyz_frac(node2,2:3)-xyz_frac(node1,2:3)).^2));
      statev_frac{i}= {frac, frac};
      N1  = CONNEC(enrElems(i),2);                                                  % Node 1 for current element
      N2  = CONNEC(enrElems(i),3);                                                  % Node 2 for current element
      N3  = CONNEC(enrElems(i),4);                                                  % Node 3 for current element
      N4  = CONNEC(enrElems(i),5);                                                  % Node 4 for current element
      X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
      Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
      xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4]; 
      [gp,gw] = gauss(2,'LINE');
      for ig = 1:length(gp)
          xi = gp(ig);                                    % Gauss points\
          statev_frac{i}{ig}.natural_coodinates = xi;
          statev_frac{i}{ig}.gauss_weight = gw(ig);
          N = 1/2*[(1-xi)   (1+xi)];
          tureCoord=[N*xynode]';
          statev_frac{i}{ig}.true_coodinates = tureCoord;
          a1 = X1+X2+X3+X4;
          a2 = -X1+X2+X3-X4;
          a3 = -X1-X2+X3+X4;
          a4 = X1-X2+X3-X4;
          b1 = Y1+Y2+Y3+Y4;
          b2 = -Y1+Y2+Y3-Y4;
          b3 = -Y1-Y2+Y3+Y4;
          b4 = Y1-Y2+Y3-Y4;
          c0 = a3*4*tureCoord(2)-a3*b1-4*tureCoord(1)*b3 + a1*b3;
          c1 = a2*b3 - a3*b2 + a4*4*tureCoord(2) - a4*b1 - 4*tureCoord(1)*b4 + a1*b4;
          c2 = a2*b4 - a4*b2;
          xi0 = roots([c2 c1 c0]);
          c0 = a2*(4*tureCoord(2) - b1)-(4*tureCoord(1) - a1)*b2;
          c1 = -a2*b3 + a3*b2 + a4*(4*tureCoord(2) - b1) - (4*tureCoord(1) - a1)*b4;
          c2 = a3*b4 - a4*b3;
          eta0 = roots([c2 c1 c0]);
          if (size(xi0,1) == 0)
              disp('Error for linear elements')
          elseif (size(xi0,1) == 1)
              xi = xi0;
          else
              ind1 = find( xi0 >= -1 );
              ind2 = find( xi0 <= 1 );
              ind = intersect(ind1,ind2);
              xi = xi0(ind); 
          end
          if (size(eta0,1) == 0)
              disp('Error for linear elements')
          elseif (size(eta0,1) == 1)
              eta = eta0;
          else
              ind1 = find( eta0 >= -1 );
              ind2 = find( eta0 <= 1 );
              ind = intersect(ind1,ind2);
              eta = eta0(ind); 
          end
          
          %double check
%           N = 1/4*[(1-xi)*(1-eta)   (1+xi)*(1-eta)  (1+xi)*(1+eta)  (1-xi)*(1+eta)];
%           compute_after = N*xyz
%           Compute_before = tureCoord
          statev_frac{i}{ig}.natural_coordinate_ele  = [xi ; eta];
      end
  end


end