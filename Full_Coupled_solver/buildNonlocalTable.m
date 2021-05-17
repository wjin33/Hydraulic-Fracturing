% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function  buildNonlocalTable(PROP)
% This function is try to build the nonlocal gauss point table 

t1=clock;

global CONNEC STATEV XYZ CRACK PSI

  nPt    = size(CRACK,1);

  NE = size(CONNEC,1);
  
  internalLength = PROP.internal_length;

  constant_c = pi*internalLength^2/3;
  
  for iElem = 1:NE
      Node = CONNEC(iElem,2:5);
      psi = PSI(Node);
      nearCrack = false;
      if (nnz(psi)==4  && abs(mean(psi)) <= internalLength)
          nearCrack = true;
      end
      xyz = XYZ(Node,2:3);
      element_center = sum(xyz)./4;
      
      [ind1,~]= find( abs(XYZ(:,2)-element_center(1)) < 1.5*internalLength ) ;
      [ind2,~]= find( abs(XYZ(:,3)-element_center(2)) < 1.5*internalLength ) ;
      ind=intersect(ind1,ind2);
      [~,~,c1] = intersect(ind,CONNEC(:,2)');                                  
      [~,~,c2] = intersect(ind,CONNEC(:,3)');                                   
      [~,~,c3] = intersect(ind,CONNEC(:,4)');                                  
      [~,~,c4] = intersect(ind,CONNEC(:,5)');                                   
      hElem    = unique([c1' c2' c3' c4']);                                       % Candidate elements for nonlocal average
      elements_group=size(hElem,2);
      for ig = 1:size(STATEV{iElem},2)
          neighbor_list=[iElem ig 1];
          local_coordintes =STATEV{iElem}{ig}.true_coodinates;
            for k = 1:elements_group
                jElem = hElem(1,k);
                for jg = 1:size(STATEV{jElem},2)
                    neighbor_coordintes=STATEV{jElem}{jg}.true_coodinates;
                    AB = neighbor_coordintes-local_coordintes;
                    length = sqrt(sum((AB).^2));
                    if length <= internalLength
                        if nearCrack
                            cross = false;
                            for ipt=1:nPt-1
                                C = CRACK(ipt,:)';
                                D = CRACK(ipt+1,:)';
                                CD = D-C;
                                AC = C - local_coordintes;
                                AD = D - local_coordintes;
                                CA = neighbor_coordintes - C;
                                CB = local_coordintes - C;
                                if (AB(1)*AC(2)-AB(2)*AC(1))*(AB(1)*AD(2)-AB(2)*AD(1)) < 0
                                    if (CD(1)*CA(2)-CD(2)*CA(1))*(CD(1)*CB(2)-CD(2)*CB(1)) < 0
                                        cross = true;
                                        break;
                                    end
                                end
                            end
                            if (~cross)
                                weight = (1- length^2/internalLength^2)^2/constant_c;
                                neighbor_list = [neighbor_list; jElem jg weight;];
                            end   
                        else
                            weight = (1- length^2/internalLength^2)^2/constant_c;
                            neighbor_list = [neighbor_list; jElem jg weight;];
                        end
                    end
                end
            end
            STATEV{iElem}{ig}.nonlocalTable = neighbor_list; 
      end
  end
  
  t2=clock;
  
  timecost = abs( t2(6) - t1(6) );
  
  X = sprintf('Build nonlocal table timecost: %f seconds', timecost);
  disp(X)

end

