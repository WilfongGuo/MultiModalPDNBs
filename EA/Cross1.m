function Offspring = Cross1(N,D,MatingPool)
Offspring=zeros(N,D);
number=sum(MatingPool,1);
P=number/size(MatingPool,1);
A=randperm(D,D);
for i=1:N
    parentsnum=randperm(size(MatingPool,1),2);
    parents1=MatingPool(parentsnum(1),:);
    parents2=MatingPool(parentsnum(2),:);
    same=parents1-parents2;
    same_p=find(same==0);
    Offspring(i,same_p)= parents1(:,same_p);
    differ=setdiff(A,same_p);
    for j=1:length(differ)
        if rand<P(differ(j))
            Offspring(i,differ(j))=0;
        else
            Offspring(i,differ(j))=1;
        end
    end

end
end
