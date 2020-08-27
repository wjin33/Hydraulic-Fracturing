function []=StateV_update()	% Initialize elastic history dependent varibles
%%
  global STATEV CONNEC NODES ELECENTER XYZ
  %
  NE = size(CONNEC, 1);
%   STATEV_new=cell(1,NE);
  %
  for iElem = 1:NE
      N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
      N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
      N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
      N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
      nodes = [N1 N2 N3 N4]';
      NN  = NODES(nodes,:);                                                   % Nodal data for current element
      ELECENTER(iElem,:) = sum(XYZ(nodes,2:3))./4;  
      HEN = nnz(NN(:,2));% Number of nodes with Heaviside enrichment
      SIGN = abs(sum(NN(:,3)));
      Ngp = size(STATEV{iElem},2);
      
      if HEN == 4 && Ngp == 4 && SIGN ~= 4
          STATEV{iElem} = GPRecovery(iElem,nodes);
%          STATEV_new{iElem} = GPRecovery(iElem,nodes);
%       else
%          STATEV_new{iElem} =  STATEV{iElem};
      end
  end
  FractureGPInitialization();
  
end

function  [Enrichedcell]=GPRecovery(iElem,nodes)
% This function recovery the statevaribles from unpropagated domain to
% propagated one

  global PSI XYZ STATEV
  
  Old_gp_coordinate = [STATEV{iElem}{1}.true_coodinates';
                       STATEV{iElem}{2}.true_coodinates';
                       STATEV{iElem}{3}.true_coodinates';
                       STATEV{iElem}{4}.true_coodinates';];
  N1 = nodes(1,1); N2 = nodes(2,1); N3 = nodes(3,1); N4 = nodes(4,1);                
  history.natural_coodinates = [0; 0];
  history.gauss_weight = 0;
  history.true_coodinates = [0; 0];
  history.sigma = [0; 0; 0; 0;];
  history.strain = [0; 0; 0; 0;];
  history.damage = [0; 0];
  history.kappa = [0; 0];
  history.EquivStrain = [0; 0];
  history.nonlocalTable = [];
  history.volume = 0;                                     % Define the global K
  
  Enrichedcell = {history, history, history, history,history, history, history, history,history, history, history, history};
  X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
  Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
  xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix
  if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
  end
  [gp,gw] = subDomain_(3,PN,xyz);                         % Full Heaviside enrichment
  
  for ig = 1:length(gp)
      xi = gp(ig,1); eta = gp(ig,2);                                    % Gauss points
      N = 1/4*[(1-xi)*(1-eta)   (1+xi)*(1-eta)  (1+xi)*(1+eta)  (1-xi)*(1+eta)];          % Shape functions
      Enrichedcell{ig}.true_coodinates = [N*xyz]';
      dist = zeros(4,1);
      for j =1:4
          dist(j,1) = sqrt(sum((Enrichedcell{ig}.true_coodinates - Old_gp_coordinate(j,:)').^2));
      end
      [~,index] = min(dist);
      Enrichedcell{ig} = STATEV{iElem}{index};
      Enrichedcell{ig}.natural_coodinates = [xi;eta];
      Enrichedcell{ig}.gauss_weight = gw(ig,1);
      Enrichedcell{ig}.true_coodinates = [N*xyz]';   
  end

end


function  FractureGPInitialization()
% This function calculates the global stiffness matrix for the desired 
% discontinuities defined by the user supplied input.

  global CONNEC NODES PSI XYZ connec_frac xyz_frac statev_frac PROP CRACK

  fracxy = [];
  enrElems = [];
  
  CRACK

  for iElem = 1:size(CONNEC,1)
      N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
      N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
      N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
      N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
      NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
      HEN = nnz(NN(:,2)); 
      SIGN = abs(sum(NN(:,3))); % Number of nodes with Heaviside enrichment
    
      if (HEN == 4) && SIGN ~= 4
          % Unenriched nodes                                                              % Enriched element
          X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
          Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
          xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix

          if numel(PSI) == 0, PN = [0 0 0 0]; else
                  PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
          end
%           if size(find( PN > 0),2)==4 || size(find( PN < 0),2)==4
%               continue;
%           end
          [~,~,~,fracnode] = subDomain(3,PN,xyz);                         % Full Heaviside enrichment

          if size( unique(fracnode,'rows'),1) ==1
              continue;
          end  
          fracxy = [fracxy; fracnode];
          enrElems = [enrElems iElem];
      end
  end
  
  frac.natural_coodinates = 0;
  frac.gauss_weight = 0;
  frac.true_coodinates = [0; 0];
  frac.natural_coordinate_ele = [0;0];
  frac.width = 0.000001;         %mm, tentatively, should be updated from mechanical calculation
  frac.fluidVelocity = 0;
  frac.seglength = 0;
  frac.jump = [0;0];
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

  UnorderedNodes = unique(fracxy,'rows');
  orderedNodes = [];
  
  k = dsearchn(UnorderedNodes,CRACK(1,:));
 % k = dsearchn(UnorderedNodes,delaunayn(UnorderedNodes),CRACK(1,:));
  orderedNodes = UnorderedNodes(k,:);
  UnorderedNodes = setdiff(UnorderedNodes,orderedNodes,'rows');
  while (size(UnorderedNodes,1)>2)
      lambda = size(orderedNodes,1);
%      k = dsearchn(UnorderedNodes,delaunayn(UnorderedNodes),orderedNodes(lambda,:));
      k = dsearchn(UnorderedNodes,orderedNodes(lambda,:));
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
  
  
%   nPt = size(CRACK,1);
%   for iSeg = 1:(nPt-1)
%       x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
%       x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
%       vectorA = [CRACK(iSeg+1,1)-CRACK(iSeg,1), CRACK(iSeg+1,2)-CRACK(iSeg,2)];
%       TempNodes = [];
%       lambda = [];
%       for inode = 1:size(UnorderedNodes,1)
%           vectorB = [UnorderedNodes(inode,1)-CRACK(iSeg,1), UnorderedNodes(inode,2)-CRACK(iSeg,2)];
%           if abs(vectorA(1)) <= 1e-6
%               if abs(vectorB(1)) <= 1e-6
%                   lambday = vectorB(2)/vectorA(2);
%                   if lambday <=1
%                     TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                     lambda = [lambda; lambday];
%                   end
%               end
%           elseif abs(vectorA(2)) <= 1e-6
%               if abs(vectorB(2)) <= 1e-6
%                   lambdax = vectorB(1)/vectorA(1);
%                   if lambdax <=1
%                     TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                     lambda = [lambda; lambdax];
%                   end
%               end
%           else     
%             lambdax = vectorB(1)/vectorA(1);
%             if lambdax <= 1
%                 if acos((vectorB(1)*vectorA(1)+vectorB(2)*vectorA(2))/(sqrt(sum(vectorB.^2))*sqrt(sum(vectorA.^2)))) <= 5
%                     TempNodes = [TempNodes; UnorderedNodes(inode,:)];
%                     lambda = [lambda; lambdax];
%                 end
%             end
%           end
%       end
%       if size(lambda,1)~=0
%           [~,I] = sort(lambda);
%           UnorderedNodes = setdiff(UnorderedNodes,TempNodes,'rows');
%           orderedNodes = [orderedNodes ; TempNodes(I,:)];
%       end
%   end
  
  segment_number = size(enrElems,2);
  xyz_frac = zeros(segment_number+1,4);
  xyz_frac(:,1) = [1:segment_number+1 ]';
  if size(orderedNodes,1) == size(xyz_frac,1)
      xyz_frac(:,2:3)  = orderedNodes;
  else
      xyz_frac(:,2:3)  = backup_UnorderedNodes;
  end
%   xyz_frac(:,2:3)  = sortrows(UnorderedNodes,2);
  
  statev_frac_new = cell(1,segment_number);
  segment_belongings = connec_frac(:,1);
  connec_frac = zeros(segment_number,3);

%   xyz_frac(:,2:3)  =  UnorderedNodes;
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
      connec_frac(i,:) = [enrElems(i) sort([node1 node2])];
      frac.seglength = sqrt(sum((xyz_frac(node2,2:3)-xyz_frac(node1,2:3)).^2));
      [index] = find(segment_belongings == enrElems(i));
      if index
          statev_frac_new{i} = statev_frac{index};
      else
        statev_frac_new{i}= {frac, frac};
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
          statev_frac_new{i}{ig}.natural_coodinates = xi;
          statev_frac_new{i}{ig}.gauss_weight = gw(ig);
          N = 1/2*[(1-xi)   (1+xi)];
          tureCoord=[N*xynode]';
          statev_frac_new{i}{ig}.true_coodinates = tureCoord;
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
          statev_frac_new{i}{ig}.natural_coordinate_ele  = [xi ; eta];
        end
      end
  end
  
  statev_frac = statev_frac_new;
  
end


function [gp,gw] = subDomain_(npt,psi,xyz)
% This function subdivides enriched elements and determines the guass 
% points and weights to be used in the integration during the assembly of 
% the stiffness matrix.
%%
corner = [1 2 3 4 1];
node   = [-1 -1;1 -1;1 1;-1 1];

% Loop through crack level set function
if isempty(psi) == 0
    for i = 1:4
        n1 = corner(i);
        n2 = corner(i+1);
        if psi(n1)*psi(n2) < 0
            r    = psi(n1)/(psi(n1)-psi(n2));
            pnt  = (1-r)*node(n1,:)+r*node(n2,:);
            xi   = pnt(1); eta = pnt(2);
            N    = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);...
                        (1+xi)*(1+eta);(1-xi)*(1+eta)];
            xpnt = dot(N,xyz(1:4,1)');
            ypnt = dot(N,xyz(1:4,2)');
            xyz  = [xyz;xpnt ypnt];
            node = [node;pnt];
        end
    end
end

% Find the triangles defining the subdomains
warning off MATLAB:delaunayn:DuplicateDataPoints
tri = delaunay(node(:,1),node(:,2));

% Loop over subtriangles to get quadrature points and weights
[q,w] = gauss(npt,'TRI');

pt = 1;
for e = 1:size(tri,1)
    coord = node(tri(e,:),:);
    xyzl  = xyz(tri(e,:),:);
    for i = 1:length(w)
        xi = q(i,1); eta = q(i,2);
        N  = [1-xi-eta;xi;eta];
        gp(pt,:) = N'*coord;
        gw(pt,1) = w(i)/2;
        pt = pt+1;
    end
end

end
