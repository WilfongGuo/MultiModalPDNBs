function [pop,tt]=FIX2(Population,functionvalue,g,D,D_score,test_adjacency)
tt= functionvalue(:,2)==0;
tt=find(tt==1);
%% ��f2=0��������޸�
for i=1:size(tt,1)
    node_position=find(Population(tt(i),:)~=0); %�ҵ���Ӧ����ѡ�����ĸ�����
    if isempty(node_position)
        Population(tt(i),:)=creatpop(1,D,D_score,g); %����ǿյ� ֱ�������µĸ���
    else
        if size(node_position,2)==1
            select=node_position;
            d=find(test_adjacency(select,:)~=0);% �ҵ���select���������Ļ���
            Can_d=D_score(d);
            Min=min(Can_d);
            index= Can_d==Min;
            add_gen=d(index');
            Population(tt(i),add_gen)=1;
        else
            %% �������㶼û���ӣ���ѡ�������������бȽ϶�����һ�������޲�
            p=randperm(size(node_position,2),2);
            node_position=node_position(:,p);
            %% �ж���������˭�ķ����󣬷�������п���Ϊ��һ��ѡ����weight��׼��
            if D_score(node_position(1))>D_score(node_position(2))
                select=node_position(1);
            else
                select=node_position(2);
            end
            %% ��ѡ�л���Ѱ�����������Ļ���
            d=find(test_adjacency(select,:)~=0);% �ҵ���select���������Ļ���
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


