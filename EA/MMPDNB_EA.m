function [PF,PS,Non_dominated_sol]=MMPDNB_EA(data,popnum,Max_CalNum,experiment_num)
%%     input:      
%                   data            :    patient sample's PGIN
%                   popnum          :    population size
%                   Max_CalNum      :    the maximum number of function evaluation
%                   experiment_num  :    the number of algorithm rus
%       output:     
%                   PF   :    the value of objective function of DNB obtained by MMPDNB
%                   PS   :    non dominated DNB obtained by MMPDNB
%                   NDS  :    the results of 30 runs of MMPDNB 

    test_adjacency=data.subnetwork_adjacency;
    D=size(test_adjacency,1); % Dimension
    test_net=zeros(D,D);
    test_net(test_adjacency~=0)=1;
    %% Calculate the degree of each gene
    Cons=CONS(test_net);
    Cons(all(Cons==0,2),:)=[];
    cons=Cons;
    Cnum=size(cons,1);
    g=cell(Cnum,1);
    for i=1:Cnum
        g{i}=find(cons(i,:)==1);
    end
    Dimension=eye(D,D);
    CV_num=Calcons(D,Cnum,g,Dimension);
    Gene_degree=abs(CV_num-Cnum);% the degree of each gene
    %% Set experiment parameters
    Non_dominated_sol=cell(experiment_num,1);
    
    for EXP_NUM=1:experiment_num
        %% Generate intial population.
        CalNum=0;   % the number of function evaluation
        Population=creatpop(popnum,D,Gene_degree,g);
        [functionvalue,calnum]=Calfunctionvalue(Population,test_adjacency); % Calculate the objective function
        CalNum=CalNum+calnum;
        [Population,tt]=FIX2(Population,functionvalue,g,D,Gene_degree,test_adjacency);
        [functionvalue,calnum]=Calfunctionvalue_afterfix(Population,test_adjacency,tt,functionvalue);
        CalNum=CalNum+calnum;
        %% Non-dominated sorting
        [FrontNo,SpCrowdDis] = M_non_domination_scd_sort(Population,functionvalue);
        
        %% Evolution
        while(CalNum<Max_CalNum)
            %% Tournament Selection
            MatingPool = TournamentSelection_hamming(popnum/2,popnum,Population,FrontNo,SpCrowdDis);
            %% Generate Offspring
            Offspring  = Cross1(popnum,D,MatingPool);
            Offspring = Muation(Offspring,Gene_degree,test_adjacency);
            [functionvalue_Offspring,calnum]=Calfunctionvalue(Offspring,test_adjacency);
            CalNum=CalNum+calnum;
            [Offspring,tt]=FIX2(Offspring,functionvalue_Offspring,g,D,Gene_degree,test_adjacency);
            [functionvalue_Offspring,calnum]=Calfunctionvalue_afterfix(Offspring,test_adjacency,tt,functionvalue_Offspring);
            CalNum=CalNum+calnum;
            %% Merge population and Offspring
            Population_new=[Population;Offspring];
            functionvalue_new=[functionvalue;functionvalue_Offspring];
            [FrontNo_new,SpCrowdDis_new] = M_non_domination_scd_sort(Population_new,functionvalue_new);
            %% Tournament Selection
            [Population,FrontNo,functionvalue,SpCrowdDis] = EnvironmentalSelection(Population_new,popnum,FrontNo_new,SpCrowdDis_new,functionvalue_new);
        end
        
        %% Record non-dominated DNB
        outputpop=Population((FrontNo==1)',:);
        Non_dominated_sol{EXP_NUM}=outputpop;
    end
    
    Non_dominated_sol=cell2mat(Non_dominated_sol);
    
    [pop,~,~]=unique(Non_dominated_sol,'rows');% De-duplication
    
    functionvalue = Calfunctionvalue(pop,test_adjacency);
    
    [FrontNo,~] = M_non_domination_scd_sort(pop,functionvalue);
    
    POP=pop(FrontNo==1,:);
    
    FV=functionvalue(FrontNo==1,:);
    
    [PF,mod_position]=sortrows(FV);
    
    PS=POP(mod_position,:);
    
    PF(:,2)=-PF(:,2);
    
end



function A_adjacent=CONS(test_Net)
[z1,z2]=find(triu(test_Net)~=0);
z=[z1,z2];
NNN=length(test_Net);

N1=NNN;
[N2,~]=size(z);
%calculate the adjacency matrix of bipartite graph
A_adjacent=zeros(N2,N1);
for i=1:N2
    
    A_adjacent(i,z(i,1))=1;
    A_adjacent(i,z(i,2))=1;
    
end
end
function CVvalue = Calcons(popnum,Cnum,g,pop)
cv=zeros(Cnum,1);
CVvalue=zeros(popnum,1);
for i=1:popnum
    ind=pop(i,:);
    for j=1:Cnum
        cv(j)=max(1-sum(ind(g{j})),0);
    end
    CVvalue(i)=sum(cv);
end
end
function pop=creatpop(popnum,D,D_score,g)
%% ѡ������һ�������Ļ���
pop=zeros(popnum,D);
gg=cell2mat(g);
%% ʣ�µ����ѡ��
for i=1:popnum
    mustselectnum_position=randperm(size(gg,1),1);  
    mustselectnum=gg(mustselectnum_position,:);
    pop(i,mustselectnum)=1;
    D_score(mustselectnum)=0;
    canditateD=find(D_score~=0);
    
    min_D_score=D_score;
    gennum=round(0.9*rand*size(canditateD,1));
    for j=1:gennum
        canditateD=find(min_D_score~=0);   % The selected gene is not participating in the selection
        variables=randperm(size(canditateD,1),2);
        if D_score(canditateD(variables(1)))>D_score(canditateD(variables(2)))
            pop(i,canditateD(variables(1)))=1;
            
            min_D_score(canditateD(variables(1)))=0;
        else
            pop(i,canditateD(variables(2)))=1;
            
            min_D_score(canditateD(variables(2)))=0;
        end
    end
end
end
function [Functionvalue,calnum] = Calfunctionvalue(pop,test_adjacency)

Functionvalue=zeros(size(pop,1),2);
functionvalue=zeros(size(pop,1),4);
for i=1:size(pop,1)
    Functionvalue(i,1)=sum(pop(i,:));
    % Calculate internal links
    a=find(pop(i,:)==1);
    matrix=test_adjacency(a,:);
    inmatrix=matrix(:,a);
    genin=nonzeros(inmatrix);
    
    functionvalue(i,1)=abs(mean(genin));
    functionvalue(i,2)=std(genin);
    
    % Calculate external links
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);
    
    functionvalue(i,3)=abs(mean(genout));
    functionvalue(i,4)=functionvalue(i,1)*functionvalue(i,2)/functionvalue(i,3);
    
end
calnum=i;
Functionvalue(:,2)= -round(functionvalue(:,4)*100)/100;

end
function [Functionvalue,calnum] = Calfunctionvalue_afterfix(pop,test_adjacency,tt,functionvalue)

for i=1:size(tt,1)
    functionvalue(tt(i),:)=Calfunctionvalue(pop(tt(i),:),test_adjacency);
end
calnum=i;
Functionvalue=functionvalue;
tt= isnan(Functionvalue(:,2));
Functionvalue(tt,2)=0;
ttt= Functionvalue(:,2)==0;
Functionvalue(ttt,1)=size(test_adjacency,2);
end

function Parents = TournamentSelection_hamming(K,N,Population,FrontNo,SpCrowdDis)
index=zeros(K,1);

hamming_dist=pdist2(Population,Population,'hamming');
[~,site2]=sort(hamming_dist,2);
K1=randperm(N,K);  % parent1 number
K2=site2(K1,2);    % parent2 number

for i=1:K
    if FrontNo(K1(i))<FrontNo(K2(i))
        index(i)=K1(i);
    elseif FrontNo(K1(i))==FrontNo(K2(i))
        if SpCrowdDis(K1(i))>=SpCrowdDis(K2(i))
            index(i)=K1(i);
        else
            index(i)=K2(i);
        end
    else
        index(i)=K2(i);
    end
    
end
Parents=Population(index,:);
end
function [epop,epop_FrontNo,functionvalue,SpCrowdDis] = EnvironmentalSelection(Population,NPE,FrontNo,SpCrowdDis,functionvalue_new)
[~,uni] = unique(Population,'rows');
SpCrowdDis     = SpCrowdDis(:,uni);
Population = Population(uni,:);
FrontNo =FrontNo(:,uni);
functionvalue_new=functionvalue_new(uni,:);
fnum=0;                                                                
while numel(FrontNo,FrontNo<=fnum+1)<=NPE                      
    if numel(FrontNo,FrontNo<=fnum+1)==numel(FrontNo,FrontNo<=fnum)
        fnum=fnum-1;
        break
    end
    fnum=fnum+1;
end
MaxFNo=fnum+1;
Next = FrontNo < MaxFNo;
Last     = find(FrontNo==MaxFNo);
[~,Rank] = sort(SpCrowdDis(Last),'descend');
Next(Last(Rank(1:NPE-sum(Next)))) = true;
epop=Population(Next,:);
epop_FrontNo=FrontNo(Next);
functionvalue=functionvalue_new(Next,:);
SpCrowdDis=SpCrowdDis(:,Next);
end
