4% close;
% clear;
% clc;
func_num=30;
func_id=1;
D=30;
Xmin=-100;
Xmax=100;
pop_size=100;
iter_max=500;
runs=25;
fhd=str2func('cec17_func');
fbest=zeros(func_num,runs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MHCHPSO
static_MHCHPSO=zeros(func_num,2);
xbest_MHCHPSO=cell(func_num,1);
goals_MHCHPSO=cell(func_num,runs);
curve_MHCHPSO=zeros(func_num,iter_max);
FEs_MHCHPSO=zeros(func_num,runs);
for i=1:func_num
    func_id=i;
    xbest=zeros(runs,D);
    curve_temp=zeros(1,iter_max);
    temp_FEs=zeros(1,runs);
    parfor j=1:runs
        % gbest     Best Solution
        % gbestval  Best Fitness
        % FEs       Function Evaluations
        % curve     The optimal fitness for each iteration (used to draw the curve) 
        [gbest,gbestval,FEs,curve]= MHCHPSO(fhd,D,pop_size,iter_max,Xmin,Xmax,func_id);
        xbest(j,:)=gbest;
        fbest(i,j)=gbestval;
        curve_temp=curve_temp+curve;
        temp_FEs(1,j)=FEs;
    end
    static_MHCHPSO(i,:)=[mean(fbest(i,:)),sum((fbest(i,:)-mean(fbest(i,:))).^2)/runs];
    xbest_MHCHPSO{i,1}=xbest;
    curve_MHCHPSO(i,:)=curve_temp/runs;
    FEs_MHCHPSO(i,:)=temp_FEs;
end
if exist(['static/',num2str(30),'D/temp'],'dir')==0
    mkdir(['static/',num2str(30),'D/temp']);
end
if exist(['xbest/',num2str(30),'D/temp'],'dir')==0
    mkdir(['xbest/',num2str(30),'D/temp']);
end
if exist(['curve/',num2str(30),'D/temp'],'dir')==0
    mkdir(['curve/',num2str(30),'D/temp']);
end
save(['static/',num2str(D),'D/temp/static_MHCHPSO.mat'],"static_MHCHPSO");
save(['xbest/',num2str(D),'D/temp/xbest_MHCHPSO.mat'],"xbest_MHCHPSO");
save(['curve/',num2str(D),'D/temp/curve_MHCHPSO.mat'],"curve_MHCHPSO");

%% Calculate the theoretical optimal solution of CEC2017 
best_fitness_cec2017=zeros(func_num,1);
for i=1:func_num
eval(['load input_data/shift_data_' num2str(i) '.txt']);
eval(['O=shift_data_' num2str(i) '(1,1:10);']);
best_fitness_cec2017(i,1)=cec17_func(O',i);
end