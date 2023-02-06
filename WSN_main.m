clear all
addpath(genpath(pwd))
func_num=1;

% test env 1
% L=100;
% R = 10;
% spacing = 1;
% dim = 2;
% point_num= 40 ;
% env_id = 1;

% test env 2
% dim=2;
% L=40;
% R=3;
% point_num=30;
% spacing = 1;
% env_id = 2;

% test env 3
dim=2;
L=80;
R=5;
point_num=30;
spacing = 1;
env_id = 3;

Xmin=0;
Xmax=L;
D=point_num*dim;
Max_FEs = 50000;
pop_size=40;
iter_max=floor((Max_FEs-1)/pop_size)+1;
runs=30;
wsn_obj = WSNC(L,R,spacing,dim,point_num);
fhd=@wsn_obj.wsn2d;
% wsn_obj.plot_panel(rand(1,80)*L);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MHCHPSO
static_MHCHPSO=zeros(func_num,2);
xbest_MHCHPSO=cell(func_num,1);
curve_MHCHPSO=zeros(func_num,iter_max);
FEs_MHCHPSO=zeros(func_num,runs);
fbest=zeros(func_num,runs);
for i=func_num:func_num
    xbest=zeros(runs,D);
    curve_temp=zeros(1,iter_max);
    temp_FEs=zeros(1,runs);
    parfor j=1:runs
        [gbest,gbestval,FEs,curve]= MHCHPSO(fhd, D, pop_size, iter_max, Xmin,Xmax, env_id);
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
fbest_MHCHPSO=fbest;
% save(['static/',num2str(D),'D/temp/static_MHCHPSO.mat'],"static_MHCHPSO");
% save(['xBest/',num2str(D),'D/temp/xbest_MHCHPSO.mat'],"xbest_MHCHPSO");
% save(['curve/',num2str(D),'D/temp/curve_MHCHPSO.mat'],"curve_MHCHPSO");
% save(['fbest/',num2str(D),'D/temp/fbest_MHCHPSO.mat'],"fbest_MHCHPSO");

%%
% static_data=zeros(6,3);
% static_data(:,2) = 1-[static_SPSO(1,1);static_CS(1,1);static_SMA(1,1);static_GNDDE(1,1);1-static_IGWO(1,1);static_MHCHPSO(1,1)];

%%
plot_coverage(wsn_obj,xbest_MHCHPSO{1},fbest_MHCHPSO,runs, env_id, "MHCHPSO")

% ------------------------¸²¸ÇÍ¼----------------------------%
function plot_coverage(wsn_obj,xbest,fbest,runs,env_id,alg_name)
    mid_val = floor((runs-1)/2)+1;
    
    [~,sort_index] = sort(fbest,1);
    wsn_obj.plot_panel(xbest(sort_index(mid_val),:), env_id, alg_name);
    
end