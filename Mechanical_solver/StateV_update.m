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
      HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
      Ngp = size(STATEV{iElem},2);
      
      if HEN == 4 && Ngp == 4 
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
  [gp,gw,~] = subDomain(3,PN,xyz);                         % Full Heaviside enrichment
  
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

  global CONNEC NODES PSI XYZ connec_frac xyz_frac statev_frac CRACK

  fracxy = [];
  enrElems = [];

  for iElem = 1:size(CONNEC,1)
      N1  = CONNEC(iElem,2);                                                  % Node 1 for current element
      N2  = CONNEC(iElem,3);                                                  % Node 2 for current element
      N3  = CONNEC(iElem,4);                                                  % Node 3 for current element
      N4  = CONNEC(iElem,5);                                                  % Node 4 for current element
      NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element             
      HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
    
      if (HEN == 4)  
          % Unenriched nodes                                                              % Enriched element
          X1 = XYZ(N1,2); X2 = XYZ(N2,2); X3 = XYZ(N3,2); X4 = XYZ(N4,2);     % Nodal x-coordinates
          Y1 = XYZ(N1,3); Y2 = XYZ(N2,3); Y3 = XYZ(N3,3); Y4 = XYZ(N4,3);     % Nodal y-coordinates
          xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                % Nodal coordinate matrix

          if numel(PSI) == 0, PN = [0 0 0 0]; else
                  PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
          end
          [~,~,~,fracnode] = subDomain(3,PN,xyz);                         % Full Heaviside enrichment
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
  segment_number = size(enrElems,2);
  statev_frac_new = cell(1,segment_number);
  old_segment_num = size(statev_frac,2);
  statev_frac_new(1,1:old_segment_num) = statev_frac;
  connec_frac = zeros(segment_number,3);
  xyz_frac = zeros(segment_number+1,4);
  maxNodeNum = max(NODES(:,2));
  xyz_frac(:,1) = [(0*maxNodeNum+1) : (0*maxNodeNum + segment_number+1)]';
  UnorderedNodes = unique(fracxy,'rows');
  orderedNodes = [];
  
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
  xyz_frac(:,2:3)  =  UnorderedNodes;
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
      frac.seglength = sqrt(sum((xyz_frac(node2,2:3)-xyz_frac(node1,2:3)).^2));
      if (i > old_segment_num)
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