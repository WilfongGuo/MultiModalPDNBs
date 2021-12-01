function [pop,tt]=FIX2(Population,functionvalue,g,D,D_score,test_adjacency)
tt= functionvalue(:,2)==0;
tt=find(tt==1);
%% 对f2=0个体进行修复
for i=1:size(tt,1)
    node_position=find(Population(tt(i),:)~=0); %找到对应个体选择了哪个基因
    if isempty(node_position)
        Population(tt(i),:)=creatpop(1,D,D_score,g); %如果是空的 直接生成新的个体
    else
        if size(node_position,2)==1
            select=node_position;
            d=find(test_adjacency(select,:)~=0);% 找到与select基因相连的基因
            Can_d=D_score(d);
            Min=min(Can_d);
            index= Can_d==Min;
            add_gen=d(index');
            Population(tt(i),add_gen)=1;
        else
            %% 即便许多点都没连接，仅选择其中两个进行比较对其中一个进行修补
            p=randperm(size(node_position,2),2);
            node_position=node_position(:,p);
            %% 判断两个基因谁的分数大，分数大的有可能为下一步选择大的weight做准备
            if D_score(node_position(1))>D_score(node_position(2))
                select=node_position(1);
            else
                select=node_position(2);
            end
            %% 对选中基因寻找与其相连的基因
            d=find(test_adjacency(select,:)~=0);% 找到与select基因相连的基因
            Can_d=D_score(d);
            Min=min(Can_d);
            index= Can_d==Min;
            add_gen=d(index');
            Population(tt(i),add_gen)=1;
        end
    end
end
pop=Population;
end


