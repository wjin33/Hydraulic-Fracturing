% Written By: Wencheng Jin, Georgia Institute of Technology (2018)
% Email: wencheng.jin@gatech.edu

function [EXTDISP,BNoset,BElset,Bsurface] = Preprocessor(file)

% Obtain the connectivity matrix and Node coordinate matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define the maximum storage space

global CONNEC NODES XYZ

maxN=50000;maxE=10000;
NODES=zeros(maxN,3);
XYZ=zeros(maxN,3);
CONNEC=zeros(maxE,5);

BNoset=cell(10,2);
BElset=cell(10,2);
Bsurface=cell(10,2);
SDISPT=[];

% fid = fopen('Nonlocal_debug.inp');
fid = fopen(file);
% fid = fopen('Arc_debug.inp');
% fid = fopen('Arc_NL_xfem.inp');
while ( ~feof(fid) )   
    tline = fgetl(fid);
    tline=strsplit(tline,',');
    if (strcmp(tline{1}, '*Node')>0)
        tline = fgetl(fid);
        tline = strsplit(tline,',');
        while(strcmp(tline{1},'*Element')<=0)
            NodeN= str2double(tline{1});
            NODES(NodeN,1)=NodeN;
            XYZ(NodeN,:)=[NodeN str2double(tline{2}) str2double(tline{3})];
            tline = fgetl(fid);
            tline = strsplit(tline,',');
        end
        XYZ=XYZ(1:NodeN,:);
        NODES=NODES(1:NodeN,:);
    end
    if (strcmp(tline{1},'*Element')>0)
        tline = fgetl(fid);
        tline = strsplit(tline,',');
        while(strcmp(tline{1},'*Nset')<=0)
            ElementN= str2double(tline{1});
            CONNEC(ElementN,:)=[ElementN str2double(tline{2}) str2double(tline{3}) str2double(tline{4}) str2double(tline{5})];
            tline = fgetl(fid);
            tline = strsplit(tline,',');
        end
        CONNEC=CONNEC(1:ElementN,:); 
    end
    if (strcmp(tline{1},'*Assembly')>0)
        tline = fgetl(fid);
        tline = strsplit(tline,',');
        Nsi=1;
        Esi=1;
        while(strcmp(tline{1},'*End Assembly')<=0)
            
            if (strcmp(tline{1},'*Nset')>0)
                [m,n]=size(tline);
                if(strcmp(tline{n},' generate')>0)
                    temp=strsplit(tline{2},'=');
                    BNoset{Nsi,1}=temp{2};
                    tline = fgetl(fid);
                    tline = strsplit(tline,',');
                    temp = [str2double(tline{1}) str2double(tline{2}) str2double(tline{3})];
                    BNoset{Nsi,2}=temp(1):temp(3):temp(2);
                    Nsi=Nsi+1;
                else
                    temp1=strsplit(tline{2},'=');
                    BNoset{Nsi,1}=temp1{2};
                    tline = fgetl(fid);
                    tline = strsplit(tline,',');
                    temp2=[];
                    while (strncmp(tline{1}, '*', 1)<=0)
                        temp=str2double(tline);
                        n=size(temp,2);
                        if isnan(temp(1,n))
                            temp=temp(1,1:n-1);
                        end
                        temp2=[temp2 temp];
                        tline = fgetl(fid);
                        tline = strsplit(tline,',');
                    end
                    BNoset{Nsi,2}=temp2;
                    Nsi=Nsi+1;
                end
            end
            
            if (strcmp(tline{1},'*Elset')>0)
                [m,n]=size(tline);
  
                if(strcmp(tline{n},' generate')>0)
                    temp1=strsplit(tline{2},'=');
                    BElset{Esi,1}=temp1{2};
                    tline = fgetl(fid);
                    tline = strsplit(tline,',');
                    temp = [str2double(tline{1}) str2double(tline{2}) str2double(tline{3})];
                    BElset{Esi,2}=temp(1):temp(3):temp(2);
                    Esi=Esi+1;
                else
                    temp1=strsplit(tline{2},'=');
                    BElset{Esi,1}=temp1{2};
                    tline = fgetl(fid);
                    tline = strsplit(tline,',');
                    while (strncmp(tline{1}, '*', 1)<=0)
                        temp=str2double(tline);
                        if (isnan(temp(1,size(temp,2))))
                            temp = temp(1);
                        end
                        BElset{Esi,2}=[BElset{Esi,2} temp];
                        tline = fgetl(fid);
                        tline = strsplit(tline,',');
                    end
                    Esi=Esi+1;
                end
            end
            if strcmp(tline{1},'*Elset')<=0 && strcmp(tline{1},'*Nset')<=0 && strcmp(tline{1},'*End Assembly')<=0
                tline = fgetl(fid);
                tline = strsplit(tline,',');
            end
        end
        Esi=Esi-1;
        Nsi=Nsi-1;
        Ssi=0;
        for i=1:Esi
            line = strsplit(BElset{i,1},'_'); 
            [~,n] = size(line);
            if n>2
                if Ssi> 0 
                    if (strcmp(Bsurface{Ssi,1},line{2})>0)
                        temp = strsplit(line{3},'S');
                        ind1 = 1:2:(2*size(BElset{i,2},2)-1);
                        ind2 = 2:2:(2*size(BElset{i,2},2));
                        temp2 = [];
                        temp2(ind2) = str2double(temp{2});
                        temp2(ind1) = BElset{i,2};
                        Bsurface{Ssi,2} = [Bsurface{Ssi,2} temp2]; 
                    else
                        Ssi = Ssi + 1;
                        Bsurface{Ssi,1}=line{2};
                        temp = strsplit(line{3},'S');
                        ind1 = 1:2:(2*size(BElset{i,2},2)-1);
                        ind2 = 2:2:(2*size(BElset{i,2},2));
                        temp2 = [];
                        temp2(ind2) = str2double(temp{2});
                        temp2(ind1) = BElset{i,2};
                        Bsurface{Ssi,2} = [Bsurface{Ssi,2} temp2];  
                    end
                else
                    Ssi = Ssi + 1;
                    Bsurface{Ssi,1}=line{2};
                    temp = strsplit(line{3},'S');
                    ind1 = 1:2:(2*size(BElset{i,2},2)-1);
                    ind2 = 2:2:(2*size(BElset{i,2},2));
                    temp2 = [];
                    temp2(ind2) = str2double(temp{2});
                    temp2(ind1) = BElset{i,2};
                    Bsurface{Ssi,2} = [Bsurface{Ssi,2} temp2];  
                end      
            end
        end
    end
    
    if (strcmp(tline{1},'*Boundary')>0)
        tline = fgetl(fid);
        tline = strsplit(tline,',');
        [m,n]=size(tline);
        while (strncmp(tline{1}, '*', 1)<=0)
            for i=1:Nsi
                if strcmp(tline{1},BNoset{i,1})>0
                    NNDB=size(BNoset{i,2},2);
                    Temp=zeros(NNDB,3);
                    Temp(:,1)=BNoset{i,2}';
                    Temp(:,2)=str2double(tline{2}).*ones(NNDB,1);
                    if n==4
                        Temp(:,3)=str2double(tline{4}).*ones(NNDB,1);
                    end
                    SDISPT=[SDISPT;Temp];   
                end
            end
            tline = fgetl(fid);
            tline = strsplit(tline,',');
        end
    end
end
fclose(fid);

TEMP = [(SDISPT(:,1)-1)*2+SDISPT(:,2)   SDISPT(:,3)];
[~,ind]=sort(TEMP(:,1));
EXTDISP = TEMP(ind,:);

end




