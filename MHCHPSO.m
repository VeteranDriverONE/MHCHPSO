function [bestX, bestFitness, FEs, curve] = MHCHPSO_FIN(fhd,dim,agent_num,max_iter,lb,ub,func_id)
%-------------------------------------------------------------------------%
%    Multi-topology hierarchical collaborative hybrid particle            %
%        swarm optimization algorithm (MHCHPSO)                           %
%                                                                         %
%    Developed in MATLAB R2018a                                           %
%                                                                         %
%    Author and programmer: Kanqi Wang                                    %
%                                                                         %
%    e-Mail: wongkq@foxmail.com                                           %
%            wongkq@stumail.nwu.edu.cn                                    %
%                                                                         %
%    Programming dates: April 2020 to August 2020                         %
%                                                                         %
%-------------------------------------------------------------------------%
%    Parameter Description :                                              %
%    Output:                                                              %
%    bestX represents the optimal solution                                %
%    bestFitness represents the optimal fitness                           %
%    FEs represents the number of function evaluations                    %
%    curve records the global optimal fitness of each iteration           %
%                                                                         %
%    Input:                                                               %
%    fhd represents the fitness function (CEC2017 is used by default )    %
%    max_iter represents the maximum number of iterations                 %
%    agent_num represents the total population number, which will         %
%        be evenly distributed to each topology                           %
%    dim represents the dimension of the problem                          %
%    lb is a vector that represents the lower bound of the search         %
%        for each dimension of the problem                                %
%    ub is a vector that represents the upper bound of the search         %
%         in each dimension of the problem                                %
%    func_id indicates the function number on CEC2017                     %              
%                                                                         %
%-------------------------------------------------------------------------%
    if ~exist('fhd','var')
        fhd=str2func('cec17_func');
        func_id=1;
        dim=30;
        agent_num=100;
        max_iter=500;
        lb=-100;
        ub=100;
    end
    if length(lb)==1
        lb=repmat(lb,1,dim);
        ub=repmat(ub,1,dim);
    end
    rng('default');
    rng('shuffle');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c=4;  % 4 topologies
    [group_best_pos, group, agent, agent_fit, velocity, history, history_fit]=init(fhd, agent_num, dim, lb, ub,func_id, c);
    t=1;
    FEs=200;
    v_bound=(ub-lb)/10;  % Maximum and minimum speed of particles 
    group_num=agent_num/c;  % Number of particles per topology 
    group_best=agent(group_best_pos,:);  % Index of the topological best solution  
    w=ones(c,1)*0.5; % Inertial weight  
    local=1; % Label of fast convergence topology 
    levy=3;  % Label of diversity topology 
    n=5; % Judgment variables of iteration difference analysis strategy (IDAS)
    state=0;
    
    group_best_fits_t=[];
    bestFitness=Inf;
    curve=zeros(1,max_iter);
    temp_state_record=zeros(max_iter,3);
    convergence_group=[];
    count_flag=1;  % Used in Fast convergence adaptive inertia weight with activation
    count=1;  % Used in Fast convergence adaptive inertia weight with activation
    p=1;  % Used in Rule3 
    while t<max_iter
        %-----------------Update collaborative topology ------------------%
        w(local,1)=rand(1)*0.6+0.2;
        w(levy,1)=0.8;
        w(2,1)=IFEntropy(group(2,:),agent,group_best_pos(2,1),v_bound,lb,ub);
        if ~isempty(convergence_group) && ismember(4,convergence_group(:,1)) % Determine whether to activate the weight 
            count_flag=t;
            count=1;
        end
        w(4,1)=0.9^(count/max_iter*50)*(-0.016*(count_flag/max_iter*50)+1);
        count=count+1;
        temp_state_record(t,4)=w(2,1);
        temp_state_record(t,3)=w(4,1);
        c1=2;
        c2=2;
        for i=1:c
            if i==local || i==levy
                continue;
            end
            velocity(group(i,:),:)=w(i,1)*velocity(group(i,:),:)...
                +c1*rand(group_num,1).*(agent(group_best_pos(i,1),:)-agent(group(i,:),:))...
                +c2*rand(group_num,1).*(history(group(i,:),:)-agent(group(i,:),:));
        end
        % -------------------Update diversity topology -------------------%
        beta=3/2;
        sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
        alpha=rand(1);
        for i=1:group_num
            s=agent(group(levy,i),:);
            u=randn(size(s))*sigma;
            v=randn(size(s));
            step=u./abs(v).^(1/beta);    
            velocity(group(levy,i),:)=alpha*step.*(agent(group_best_pos(levy,1),:)-s).*randn(size(s));
            agent(group(levy,i),:)=Bounds(history(group(levy,i),:)+velocity(group(levy,i),:),lb,ub);
        end
        %------------------Update fast convergence topology---------------%
        if state==1
             corre=lamarcking(agent_fit_old(group(local,:),1),agent_fit(group(local,:),1),corre,eta,setZero);
        else
            eta=zeros(group_num,1);
            [corre]=competition(group_num,history_fit(group(local,:)),0);
            state=1;
        end
        c1=2;
        c2=2;
        velocity(group(local,:),:)=w(local,1).*velocity(group(local,:),:)...
            +c1*rand(group_num,1).*(agent(group_best_pos(local,1),:)-agent(group(local,:),:))...
                +c2*rand(group_num,1).*(history(group(local,corre),:)-agent(group(local,:),:));
        velocity=Bounds(velocity,-v_bound,v_bound);
        for i=1:c
            if i==levy
                continue;
            end
            agent(group(i,:),:)=Bounds(agent(group(i,:),:)+velocity(group(i,:),:),lb,ub);
        end
        agent_fit_old=agent_fit;
        setZero=[];
        for i=1:agent_num
            agent_fit(i,1)=fhd(agent(i,:)',func_id);
            if agent_fit(i,1)<history_fit(i,1)
                history(i,:)=agent(i,:);
                history_fit(i,1)=agent_fit(i,1);
                if sum(find(group(local,:))==i)>0
                    setZero=[setZero;find(group(local,:)==i)];
                end
            end
            if agent_fit(i,1)<bestFitness
                bestFitness=agent_fit(i,1);
                bestX=agent(i,:);
            end
        end
        FEs=FEs+agent_num;
        %----------------Update the best of each topology ----------------%
        for i=1:c
            [~,group_min_pos]=min(agent_fit(group(i,:),1));
            if agent_fit(group(i,group_min_pos),1)<fhd(group_best(i,:)',func_id)
                %  This part can reduce the function evaluation of this part by 
                %    adding variables, so the calculation of FEs is not increased 
                group_best_pos(i,1)=group(i,group_min_pos);
                group_best(i,:)=agent(group_best_pos(i,1),:);
            else
                agent(group_best_pos(i,1),:)=group_best(i,:);
                agent_fit(group_best_pos(i,1),1)=fhd(group_best(i,:)',func_id);
                %  This part can reduce the function evaluation of this part 
                %    by adding variables, so the calculation of FEs is not increased 
            end
        end
        if size(group_best_fits_t,2)<n % Save the optimal fitness of all topologies in the last n iterations 
            group_best_fits_t=[group_best_fits_t,agent_fit(group_best_pos(:,1),1)];
        else
            group_best_fits_t=[group_best_fits_t(:,2:end),agent_fit(group_best_pos(:,1),1)];
        end
        %--------------------Rule2---------------------%
        if group_best_fits_t(local,1)-group_best_fits_t(local,end)<0.001 ...
                && length(group_best_fits_t(local,:))==n  
            % If there is no change in the optimal solution for n iterations, it is considered to have converged
            [~,temp_pos]=min(agent_fit(group_best_pos,1));
            temp_state_record(t,2)=1;
            ri=randi(group_num);
            if temp_pos~=local  % Execute Rule2
                agent(group_best_pos(local,1),:)=agent(group_best_pos(temp_pos,1),:);
                agent_fit(group_best_pos(local,1),1)=agent_fit(group_best_pos(temp_pos,1),1);
                history(group_best_pos(local,1),:)=history(group_best_pos(temp_pos,1),:);
                group_best(local,:)=agent(group_best_pos(local,1),:);
            end
            group_best_fits_t=group_best_fits_t(:,end);
        end
        %--------------------Topological interaction ---------------------%
        convergence_group=control_strategy(agent,group,group_best_pos,v_bound,[local,levy],group_best_fits_t);
        if isempty(convergence_group) 
        %  All collaborative topologies have not converged, and no operation is performed 
            temp_state_record(t,1)=0;
            excludes=1:1:c;
        elseif size(convergence_group,1)==c-length([levy,local]) 
            %  All collaborative topologies converge, choose topology to execute Rule3
            excludes=[local,levy];
            if rand(1)<p
                p=0;
                temp_state_record(t,1)=3;
                levy_member=randi(group_num,[1,c]);
                for i=1:c
                    if ismember(i,convergence_group(:,1))
                        middle_member=randi(group_num);
                        agent(middle_member,:)=agent(levy_member(1,i),:);
                        history(middle_member,:)=history(levy_member(1,i),:);
                    end
                end
            end
        else %  Any cooperative topology convergence, 
            excludes=setdiff(1:1:c,convergence_group(:,1)'); % Choose topology to execute Rule4
            if sum(convergence_group(:,2))==1 % Execute Rule1
                temp_state_record(t,1)=2;
                local_pos=group(local,randi(group_num));
                middle_best_pos=group_best_pos(convergence_group(convergence_group(:,2)==1),1);
                agent(local_pos,:)=agent(middle_best_pos,:);
                history(local_pos,:)=history(middle_best_pos,:);
                velocity(local_pos,:)=history(middle_best_pos,:);
                eta(local_pos,1)=max(eta);
            elseif sum(convergence_group(:,2))==0
                temp_state_record(t,1)=1;
            else
                disp("other situation in convergence_group");
            end
        end
        if p<1
            p=p+0.1;
        end
        index=zeros(c,1);
        for i=1:c  
            index(i,1)=randi(group_num);
            while group(i,index(i,1))==group_best_pos(i,1)
                index(i,1)=randi(group_num);
            end
        end
        for i=1:c-1 % Exchange information between topologies
            if ~ismember(i,excludes)
                continue;
            end
            cur=i;
            next=i+1;
            while ismember(next,[local,levy])
                next=next+1;
                if next>c
                    next=1;
                end
            end
            temp1=group(cur,index(cur,1));
            group(cur,index(cur,1))=group(next,index(next,1));
            group(next,index(next,1))=temp1;
        end
        curve(1,t)=bestFitness;
        t=t+1;
    end
end
%% Auxiliary function 
%--------------------------Boundary function -----------------------------%
function [des]=Bounds(source,lb,ub)
    flagUb=source>ub;
    flagLb=source<lb;
    des=(source.*(~(flagUb+flagLb)))+ub.*flagUb+lb.*flagLb;
end
%--------------------------Initialization function -----------------------%
function [group_best_pos, group, agent, agent_fit, velocity, history, history_fit]=init(fhd,agent_num,dim,lb,ub,func_id,c)
    agent=rand(agent_num,dim).*(ub-lb)+lb;
    history=rand(agent_num,dim).*(ub-lb)+lb;
    velocity=rand(agent_num,dim).*(ub-lb)/10-(ub-lb)/20;
    agent_fit=zeros(agent_num,1);
    history_fit=zeros(agent_num,1);
    for i=1:agent_num
        agent_fit(i,1)=fhd(agent(i,:)',func_id);
        history_fit(i,1)=fhd(history(i,:)',func_id);
    end
    % Divide the population into topologies and obtain the optimal solutions for each topology 
    [group_best_pos,group]=selectGroupBest(agent, agent_fit, c);
end
%------------------------Construct topology ------------------------------%
function [group_best_pos,group]=selectGroupBest(agent, agent_fit, c)
    % group_best_pos records the optimal solution for each topology 
    % group records each topology member 
    agent_num=size(agent,1);
    index=crossvalind('Kfold', agent_num, c); 
    group=zeros(c,agent_num/c);
    group_best_pos=zeros(c,1);
    for i=1:c
        group(i,:)=find(index==i)';
        [~,sort_index]=sort(agent_fit(group(i,:),1));
        group_best_pos(i,1)=group(i,sort_index(1,1));
    end
end
%--------------------------Lamarck mechanism -----------------------------%
function [corre]=lamarcking(agent_fit_old,agent_fit,corre,eta,setZero)
    % The value j on the i-th dimension of corre indicates that particle i 
    %     chooses the historical optimal solution of particle j as the learning goal  
    num=length(agent_fit_old);
    temp_eta=zeros(num,2);
    diff=agent_fit_old-agent_fit;
    for i=1:num
        temp_eta(corre(i,1),1)=temp_eta(corre(i,1),1)+diff(i,1);
        if diff(i,1)>0
            temp_eta(corre(i,1),2)=temp_eta(corre(i,1),2)+1;
        elseif diff(i,1)<0
            temp_eta(corre(i,1),2)=temp_eta(corre(i,1),2)-1;
        end
    end
    [~,compre_index]=sortrows([temp_eta(:,2),temp_eta(:,1)],[1,2],'ascend');
    for i=1:num
        eta(i,1)=eta(i,1)+tanh(find(compre_index==i)/num*8-4);
    end
    if ~isempty(setZero)
        eta(setZero)= eta(setZero)+0.1;
    end
    eta=max(eta,-4);
    eta=min(eta,4);
    positive_num=0;
    for i =1:num
        if eta(i,1)>0
            positive_num=positive_num+1;
            positive_eta(positive_num,:)=[eta(i,1),i];
        end
    end
    positive_eta(1:positive_num,:)';
    if positive_num<2  % Start league mechanism based on existing eta 
        [corre]=competition(num,eta,1);
    else
        [corre]=competition(num,positive_eta(1:positive_num,:),1);
        corre=positive_eta(corre,2);
    end
end
%----------------------Information exchange strategy ---------------------%
function [result]=control_strategy(agent,group,group_best_pos,v_bound,exclude,group_fits_t)
    % result record the number of collaborative topology convergence, 
    %      topology number and whether there is a full topology optimal solution 
    % Determine which interaction strategy to implement 
    % 0 means all topologies are not interactive 
    % 1 means that any cooperative population converges, 
    %      and the convergence result is not fully topologically optimal 
    % 2 means that any cooperative topology converges, 
    %      and the convergence result is the optimal solution of the whole topology 
    % 3 means that both intermediate populations have converged 
    result=[];
    [c,group_num]=size(group);
    [~,min_c]=min(group_fits_t(:,end));
    count=0;
    for i=1:c
        if ~ismember(i,exclude) % Topology i is a collaborative topology 
            dist=abs(agent(group_best_pos(i,1),:)-agent(group(i,:),:));
            num=prod(dist<v_bound,2); 
            if sum(num)>=group_num*1 && sum((group_fits_t(i,1)-group_fits_t(i,end))<0.001)==0
                % If all the topologies are in the neighborhood of the optimal solution, 
                %     and the optimal solution is changed for n iterations 
                count=count+1;
                if min_c==i
                    result(count,:)=[i,1];
                else
                    result(count,:)=[i,0];
                end
            end
        end
    end
end
%--------------------------League mechanism ------------------------------%
function [corre]=competition(competitor_num,score,sign)
    % competitor_num indicates the number of competitors 
    % score represents the value of ¦Ç 
    % sign=0, indicates that the value is small first; 
    % sign=1, indicates that the value is large first
    corre=zeros(competitor_num,1);
    num=size(score,1);
    for i=1:competitor_num
       ran_num=floor(rand(1,2)*num)+1;
       if (score(ran_num(1,1),1)<score(ran_num(1,2),1) && sign==0) || (score(ran_num(1,1),1)>score(ran_num(1,2),1) && sign==1)
           corre(i,1)=ran_num(1,1);
       else
           corre(i,1)=ran_num(1,2);
       end
    end
end
function output=sigmoid(x)
    output=1./(1+exp(-x));
end
%--------------------Intuitionistic fuzzy entropy ------------------------%
function w=IFEntropy(group,agent,best_pos,v,lb,ub)
    % w indicates inertia weight 
    group_num=size(group,2);
    dist=agent(best_pos,:)-agent(group,:);
    cluster_num=sum(prod(dist<v,2));
    miu=cluster_num/group_num;
    gamma=1-miu;
    pi=max(sqrt(sum(dist.^2,2)))/sqrt(sum((ub-lb).^2));
    miu=miu*(1-pi);
    gamma=gamma*(1-pi);
    entorpy=(min(miu,gamma)+pi)/(max(miu,gamma)+pi);
    w=0.4*rand(1)+0.8*entorpy;
end