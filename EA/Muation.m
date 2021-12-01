function Offspring = Muation(pop,D_score,test_adjacency)
Offspring=zeros(2*size(pop,1),size(pop,2));

weight=sum(test_adjacency,2);
canditate=find(D_score==1);
for i=1:size(pop,1)
    offspring=pop(i,:);
    selected=find(offspring(1,:)==1);
    
    %% 计算该解的内外矩阵
    matrix=test_adjacency(selected,:);
    inmatrix=matrix(:,selected);
    genin=nonzeros(inmatrix);
    
    gen=nonzeros(matrix);
    genout=setdiff(gen,genin);
    before=sum(genout);
    
    if rand<0.5                 %添加
        canditate_node=setdiff(canditate,selected);
        
        if isempty(canditate_node)
            Offspring((2*i)-1,:)=pop(i,:);
        else
            after=before+weight(canditate_node,:);
            after=abs(after);
            value=min(after);
            p=find(after==value);
            
            P=canditate_node(p(1));
            offspring(1,P)=1;
            Offspring((2*i)-1,:)=offspring;
            if length(p)~=1
                offspring=pop(i,:);
                P=canditate_node(p(2));
                offspring(1,P)=1;
                Offspring(2*i,:)=offspring;
            end
        end
    else                        %减少
        alone=sum(inmatrix,1);
        alone= alone==0;
        canditate_node=selected(alone);
        if isempty(canditate_node)
            Offspring((2*i)-1,:)=pop(i,:);
        else
            after=before-weight(canditate_node,:);
            after=abs(after);
            value=min(after);
            p=find(after==value);
            P=canditate_node(p(1));
            offspring(1,P)=0;
            Offspring((2*i)-1,:)=offspring;
            if length(p)~=1
                offspring=pop(i,:);
                P=canditate_node(p(2));
                offspring(1,P)=0;
                Offspring(2*i,:)=offspring;
            end
        end
    end
end
value=sum(Offspring,2);
nan=value==0;
Offspring(nan,:)=[];
end


