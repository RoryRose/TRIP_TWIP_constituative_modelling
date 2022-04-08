%% open file
close all
clear all
[file, path] = uigetfile('.xlsx','Select Material Data sheet','MultiSelect','off');
[tensfile,tenspath]=uigetfile('.xls','Select Tensile Data File','MultiSelect','off');
%% extract datafile
%read xls file as a table for easy proscessig

consts=readtable([path,file],"Sheet",1);
exp_dat=readtable([path,file],"Sheet",2);
stess_dat=readtable([tenspath,tensfile],"Sheet",1);
%% fit f_alpha
%define how we are treating the ranges in composition
%define experimental strain
e=[exp_dat.FC_Tstrain;exp_dat.FC_Tstrain];
modelstrain=[0:de:max(stess_dat.Tstrain)];
f_mart=[exp_dat.Martensite_lower;...
    exp_dat.Martensite_upper]./100;
f_austen=[exp_dat.Austenite_lower;...
    exp_dat.Austenite_upper]./100;
%% basic 
[sigma,p,f,dei,fitconst]=f_constituativemodel(consts,exp_dat,plotq,verbose,f_mart,f_austen,e,modelstrain);
%% sensitivity
de=1e-2;
verbose=0;
plotq=0;
vars=consts.Constant;

pcttotest=[5];
con_lst=combvec([1:length(vars)],[1:3],[1:length(pcttotest)]);
gof=NaN(1,length(con_lst));
sigmaT=cell(1,length(con_lst));
for idx=1:length(con_lst)
    varidx=con_lst(1,idx);
    elemidx=con_lst(2,idx);
    pctidx=con_lst(3,idx);
    consts2=consts;
    consts2(varidx,elemidx+2)=...
        {table2array(consts(varidx,elemidx+2))+...
        table2array(consts(varidx,elemidx+2))./...
        100.*pcttotest(pctidx)};
    sigmaT{idx}=f_constituativemodel(consts2,exp_dat,plotq,verbose,f_mart,f_austen,e,de);
    goft=sigmaT{idx}(end-1);
    gof(idx)=goft;%replace later with the actual gof to experimental data
end
[sigma,p,f,dei,fitconst]=f_constituativemodel(consts,exp_dat,plotq,verbose,f_mart,f_austen,e,de);


%%
pltgof=(gof(con_lst(3,:)==1)-sigma(end-1))./sigma(end-1).*100;
pltgof=reshape(pltgof,[length(vars),3]);
figure(1)
b=bar3(abs(pltgof),'detached');
set(b,'FaceAlpha',0.5)
zlabel('% change final stress fo4 5% vartiable pertubation')
ylabel('Variable')
xlabel('Element')
% Change the x and y axis tick labels
yticks([1:length(vars)])
set(gca, 'YTickLabel', vars)
set(gca, 'XTickLabel', {'Austenite','Ferrite','Martensite'})
