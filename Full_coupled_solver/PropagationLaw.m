function [Propagation] = PropagationLaw

    global CRACK ELECENTER PROP XYZ
    
    Xmin=min(XYZ(:,2));
    Xmax=max(XYZ(:,2));
    Ymin=min(XYZ(:,3));
    Ymax=max(XYZ(:,3));
    nPt    = size(CRACK,1);           % Number of data points defining crack
    
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
    
    Propagation = false ;
    
    critical_damage = 0.1;
    
    internal_length = PROP.internal_length;
    
    if nCT == 1
    
        Tip = CRACK(nPt, :);
        crack_dirction = CRACK(nPt, :) - CRACK(nPt-1, :);
    
        [I1,~] = find(abs(ELECENTER(:,1)-Tip(1)) < 1.5*internal_length );
        [I2,~] = find(abs(ELECENTER(:,2)-Tip(2)) < 1.5*internal_length );
        Tip_neighbor = intersect(I1,I2);
    
        [scale,weighted_vector,weighted_damage] = Weigt_value_in_tip(crack_dirction,Tip_neighbor,Tip,internal_length);
        
        weighted_damage = weighted_damage./scale;
        [I,~]=find(weighted_damage >= critical_damage);
        if size(I)>1
            [I,~]=find(weighted_damage == max(weighted_damage));
        end
        if (~isempty(I))
            potetial_direction_uni = weighted_vector(I,:)./ sqrt(sum(weighted_vector(I,:).^2));
            nPt = nPt + 1;
            CRACK(nPt,:) = Tip + potetial_direction_uni.*internal_length/2;
            Propagation = true ;
%           nPt = nPt + 1;
%           CRACK(nPt,:) = Tip + [1,0.0]*internal_length/5;
%           Propagation = true ;
        end
    else
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% crack end tip
        Tip = CRACK(nPt, :);
        crack_dirction = CRACK(nPt, :) - CRACK(nPt-1, :);
        [I1,~] = find(abs(ELECENTER(:,1)-Tip(1)) < 1.5*internal_length );
        [I2,~] = find(abs(ELECENTER(:,2)-Tip(2)) < 1.5*internal_length );
        Tip_neighbor = intersect(I1,I2);
        [scale,weighted_vector,weighted_damage] = Weigt_value_in_tip(crack_dirction,Tip_neighbor,Tip,internal_length);
        weighted_damage = weighted_damage./scale;
        [I,~]=find(weighted_damage >= critical_damage);
        if size(I,1)>1
            [I,~]=find(weighted_damage == max(weighted_damage));
        end
        if (~isempty(I))
            potetial_direction_uni = weighted_vector(I,:)./ sqrt(sum(weighted_vector(I,:).^2));
            nPt = nPt + 1;
            CRACK(nPt,:) = Tip + potetial_direction_uni.*internal_length/2;
            Propagation = true ;
%           nPt = nPt + 1;
%           CRACK(nPt,:) = Tip + [0,1]*internal_length/2;
%           Propagation = true ;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% crack start tip
        Tip = CRACK(1, :);
        crack_dirction = CRACK(1, :) - CRACK(2, :);
        [I1,~] = find(abs(ELECENTER(:,1)-Tip(1)) < 1.5*internal_length );
        [I2,~] = find(abs(ELECENTER(:,2)-Tip(2)) < 1.5*internal_length );
        Tip_neighbor = intersect(I1,I2);
    
        [scale,weighted_vector,weighted_damage] = Weigt_value_in_tip(crack_dirction,Tip_neighbor,Tip,internal_length);
        
        weighted_damage = weighted_damage./scale;
        [I,~]=find(weighted_damage >= critical_damage);
        if size(I)>1
            [I,~]=find(weighted_damage == max(weighted_damage));
        end
        if (~isempty(I))
            potetial_direction_uni = weighted_vector(I,:)./ sqrt(sum(weighted_vector(I,:).^2));
            nPt = nPt + 1;
            CRACK(2:nPt,:) = CRACK(1:nPt-1,:);
            CRACK(1,:) = Tip + potetial_direction_uni.*internal_length/2;
            Propagation = true ;
%             nPt = nPt + 1;
%             CRACK(2:nPt,:) = CRACK(1:nPt-1,:);
%             CRACK(1,:) = Tip + [0,-1].*internal_length/2;
%             Propagation = true ;
        end
        
    end
        
end
    
function [scale,weighted_vector,weighted_damage] = Weigt_value_in_tip(crack_dirction,Tip_neighbor,Tip,internal_length)
    
    global STATEV

    number = size (Tip_neighbor,1) ;
    
    scale = 0;
    weighted_vector = [0 0;0 0];
    weighted_damage = [0; 0];
    
    for i = 1:number
        iElem = Tip_neighbor(i);
        ngp =  size(STATEV{iElem},2);
        for j = 1:ngp
            TipToGP = STATEV{iElem}{j}.true_coodinates' - Tip;
            if ( crack_dirction*TipToGP' > 0.0 )
                r = sqrt(sum(TipToGP.^2));
                vector = TipToGP./r;
                if ( r < internal_length )
%                     weight = 1/((2*pi)^(3/2)*internal_length^3)*exp(-r^2/2/internal_length^2);   % Gauss Weight function
                    weight = ( 1 - r^2/internal_length^2)^2;   % Bellshape function to be consistant with nonlocal regularization
                    volume = STATEV{iElem}{j}.volume;
                    damage = (weight*volume).*STATEV{iElem}{j}.damage ;
                    weighted_damage = weighted_damage + damage;
                    weighted_vector = weighted_vector + [damage(1,1).*vector ;  damage(2,1).*vector];
                    scale = scale + weight*volume;
                end
            end
        end 
    end
    
end
