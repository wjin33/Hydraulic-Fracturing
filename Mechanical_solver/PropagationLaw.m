function [Propagation] = PropagationLaw

    global CRACK ELECENTER PROP STATEV
    
    Propagation = false ;
    
    critical_damage = 0.1;
    
    internal_length = PROP.internal_length;

    nPt    = size(CRACK,1);           % Number of data points defining crack
    Tip = CRACK(nPt, :);
    crack_dirction = CRACK(nPt, :) - CRACK(nPt-1, :);
    
    [I1,~] = find(abs(ELECENTER(:,1)-Tip(1)) < 1.5*internal_length );
    [I2,~] = find(abs(ELECENTER(:,2)-Tip(2)) < 1.5*internal_length );
    Tip_neighbor = intersect(I1,I2); 
    
    number = size (Tip_neighbor,1) ;
    
    scale = 0;
    weighted_vector = [0 0;];
    weighted_damage = 0;
    
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
                    weighted_vector = weighted_vector + damage.*vector;
                    scale = scale + weight*volume;
                end
            end
        end 
    end
    weighted_damage = weighted_damage./scale;
%     [I,~]=find(weighted_damage >= critical_damage);
%     if size(I)>1
%         [I,~]=find(weighted_damage == max(weighted_damage));
%     end
%     if (~isempty(I))
% %         potetial_direction_uni = weighted_vector(I,:)./ sqrt(sum(weighted_vector(I,:).^2));
% %         nPt = nPt + 1;
% %         CRACK(nPt,:) = Tip + potetial_direction_uni.*internal_length/2;
% %         Propagation = true ;
    if (weighted_damage>=critical_damage)
        nPt = nPt + 1;
        CRACK(nPt,:) = Tip + [0,1.0]*internal_length/4;
        Propagation = true ;
    end 

end