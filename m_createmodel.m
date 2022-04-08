%% open file
close all
clear all

[file, path] = uigetfile('.xlsx','Select Material Data sheet','MultiSelect','off');
[tensfile,tenspath]=uigetfile('.xls','Select Tensile Data File','MultiSelect','off');
sampname=strcat(file(1:2),' WITH Martensite Hardening');
%% extract datafile
%read xls file as a table for easy proscessig
consts=readtable([path,file],"Sheet",1);
exp_dat=readtable([path,file],"Sheet",2);
stress_dat=readtable([tenspath,tensfile],"Sheet",1);
% select range for experimental data - you must select zero point for
%plastic strain and endpoint where fracture starts
h=figure(100);
ax=gca;
plot(stress_dat.Tstrain,stress_dat.Tstress)
roi=drawline(ax);
input('Happy to procede? Type any key to continue:','s');
%define the line you just drew
pos1=roi.Position;
%define the ROI line equation and sample 1000 points on it
startstrain=min([pos1(2,1),pos1(1,1)]);
endstrain=max([pos1(2,1),pos1(1,1)]);
idx=find(stress_dat.Tstrain>=startstrain & stress_dat.Tstrain<=endstrain);
plot(stress_dat.Tstrain(idx),stress_dat.Tstress(idx))
stress_dat=stress_dat(idx,:);
stress_dat.Pstrain=stress_dat.Tstrain-stress_dat.Tstrain(1);
exp_dat.Pstrain=exp_dat.FC_Tstrain-startstrain;
figure(100)
plot(stress_dat.Pstrain,stress_dat.Tstress)
pause(0.2)
close(100)
%% define how we are treating the ranges in composition and range for plastic experimental data
%=====================================
%which datapoints to use for volume fractions?
e=[exp_dat.Pstrain;exp_dat.Pstrain];
f_mart=[exp_dat.Martensite_lower;...
    exp_dat.Martensite_upper]./100;
f_austen=[exp_dat.Austenite_lower;...
    exp_dat.Austenite_upper]./100;
% alternative for debugging
%{
f_mart=exp_dat.Martensite_upper./100;
f_austen=exp_dat.Austenite_lower./100;
e=exp_dat.FC_Tstrain;
%}
%=================================================
%User defined variables

%consts.martensite(strcmp(consts.Constant,'k1'))=0;
%consts.martensite(strcmp(consts.Constant,'k2'))=0;

% model f_alpha using base values
close all
%==================================================
%settings structure for the model

settings=struct();
settings.plotq=1;
settings.de=1e-2;
settings.verbose=1;
settings.vars=consts.Constant;
%decide what you want and how to model each component

%MODEL COMPONENTS TO INCLUDE
settings.twipmodel = 1;%include twin spacing in dislocation mean free path? 1 for yes
settings.tripmodel = 1;%include martensite lath size in dislocation mean free path? 1 for yes
settings.austenmodel=1;%include austenite grain size in dislocation mean free path? 1 for yes
settings.ferritemodel=1;%include ferrite grain size in dislocation mean free path? 1 for yes
settings.ferrite = 1;%include ferrite stress contribution? 1 for yes
settings.austen = 1;%include austenite stress contribution? 1 for yes
settings.marten = 1;%include martensite stress contribution? 1 for yes


%FITTING EXPERIMENTAL PHASE FRACTIONS
settings.f_gamma_opt = 4;%1- do linear interp so use the individual sample points, 
% 2 - do linear fit of all data, 3 - just use the first value, 4 - take
% remainder from austenite and martensite, 5 - also take away twin fraction
settings.f_a_opt = 2;% 1 - take remainder form other components, 2 - just use first value
settings.f_a2_opt = 3;% 1- follow paper fitting for alpha and beta, 2 - Linear Interpolant, 3 - Avrami Model
settings.f_T_opt=2;%1 - f_austenite evaluated for first point, 2 - f_austenite evaluated at every point
%This also changes the twin model which will automaticaly assume the deafults below:
settings.deafult_alpha = 3.5;%deafult value of alpha - CHECK THIS
settings.deafult_beta = 1.1;%deafult value of beta - CHECK THIS
%=================================================

%Define user data structure

% create structure of strains
modelstrain=0:settings.de:max(stress_dat.Pstrain);
user_dat=struct();
user_dat.modelstrain=modelstrain;
user_dat.e=e-startstrain;
user_dat.f_mart=f_mart;
user_dat.f_austen=f_austen;
%==============================================
%% Run the model
[sigma,p,f,dei,~,compS]=f_constituativemodel(consts,exp_dat,user_dat,settings);
% EXPERIMENTAL DATA
p_comparison(modelstrain,compS,f,sigma,stress_dat,user_dat,consts,settings,exp_dat)
%% Try a range of f0 to see how the fit 
%{
%=======================================
%SETTINGS for fit
resetq=0;
adjustvar='a';%VARIABLE TO FIT
phase='martensite';%PHASE OF VARIABLE TO FIT
pctvariation=20;%PERCENT VARIATION IN VALUE (WILL BE SYMMETRIC ABOUT THE CURRENT VALUE)
numvals=10;%Number of values to try
user_dat.Sunc=1e6;%estimated uncertainty in stress (Pa)
user_dat.Eunc=1e-6;%estimated uncertainty in strain
settings.fitq=2;%What do you want to fit to? 1 - stress strain curve, 2 - strain hardening curve
%=====================================
%Create anonymous Function
settings3=settings;
settings3.plotq=2;
settings3.verbose=-1;
if resetq == 1
    consts2 = consts;
else
    if exist('constsFIT',"var")
        consts2=constsFIT;%ONLY USE IF DOING MULTIPLE OPTIMISATIONS!
    else
       consts2 = consts;
    end
end
startpoint=consts.(phase)(strcmp(consts.Constant,adjustvar));
%{
%ONLY FOR CHANGING A OR B
startpoint=11;
%consts.austenite(strcmp(consts.Constant,'b_prime'))=3;
%}

fun=@(var)fit_constituative_model(var,consts2,exp_dat,user_dat,settings3,stress_dat,adjustvar,phase);

%{

%}
valstotest=linspace(startpoint-(startpoint./100.*pctvariation),startpoint+(startpoint./100.*pctvariation),numvals);
%{
%set up the fit

%fit to the data
opt=optimoptions("fmincon","FiniteDifferenceStepSize",1e-3,"MaxFunctionEvaluations",20,"Display","iter-detailed")
lb = 0;
ub = 1;
f0_fit=fmincon(f_f0,startpoint,[],[],[],[],lb,ub);
%}
% OR

gof=NaN(1,length(valstotest));
for idx=1:length(valstotest)
    disp(strcat('Iteration',{' '},num2str(idx)))
    gof(idx)=fun(valstotest(idx));
end
idmin=find(gof==min(gof));
fo_fit=valstotest(idmin(1));
figure(4)
plot(valstotest,gof)
xlabel(strcat(adjustvar,{' ( '},consts.Units(strcmp(consts.Constant,adjustvar)),' )'));
ylabel('\chi^2')
settings2=settings;
settings2.plotq=2;
fit_constituative_model(fo_fit,consts2,exp_dat,user_dat,settings2,stress_dat,adjustvar,phase);
constsFIT=consts2;
constsFIT.(phase)(strcmp(consts2.Constant,adjustvar))=fo_fit;
%}
%% fit multiple variables using fminunc
%=======================================
%SETTINGS for fit
resetq=1;
adjustvar={'a','b_prime'};%VARIABLE TO FIT
phase={'martensite','martensite'};%PHASE OF VARIABLE TO FIT
numvals=100;%Number of values to try
user_dat.Sunc=1e6;%estimated uncertainty in stress (Pa)
user_dat.Eunc=1e-3;%estimated uncertainty in strain
settings.fitq=2;%What do you want to fit to? 1 - stress strain curve, 2 - strain hardening curve
%=====================================
%Create anonymous Function
settings3=settings;
settings3.plotq=0;
settings3.verbose=-1;
consts.martensite(strcmp(consts.Constant,'b_prime'))=3;
consts.martensite(strcmp(consts.Constant,'a'))=3;
if resetq == 1
    consts2 = consts;
else
    if exist('constsFIT',"var")
        consts2=constsFIT;%ONLY USE IF DOING MULTIPLE OPTIMISATIONS!
    else
       consts2 = consts;
    end
end
startpoint=NaN(1,length(adjustvar));
for i=1:length(adjustvar)
    startpoint(i)=consts.(phase{i})(strcmp(consts.Constant,adjustvar{i}));
end

fun=@(var)fit_constituative_model(var,consts2,exp_dat,user_dat,settings3,stress_dat,adjustvar,phase);

%fit to the data
opt=optimoptions("fmincon","MaxFunctionEvaluations",numvals,"PlotFcn","optimplotfval");
f0_fit=fmincon(fun,startpoint,[],[],[],[],[0,1],[inf,inf],[],opt);
settings2=settings3;
settings2.plotq=2;
fit_constituative_model(f0_fit,consts2,exp_dat,user_dat,settings2,stress_dat,adjustvar,phase);
constsFIT=consts2;
for i=1:length(adjustvar)
    constsFIT.(phase{i})(strcmp(consts.Constant,adjustvar{i}))=f0_fit(i);
end
%% save the data in a excel sheet and save the figures
filename=strcat('Summary Data from fitted ',sampname,' data WITH twinning');
settings2=settings;
settings2.twipmodel=1;
settings2.plotq=0;
[sigma,p,f,dei,~,compS]=f_constituativemodel(constsFIT,exp_dat,user_dat,settings2);

p_comparison(modelstrain,compS,f,sigma,stress_dat,user_dat,constsFIT,settings2,exp_dat)
h=figure(3);
saveas(h,filename,'tiffn')
saveas(h,filename,'fig')
etot=[user_dat.modelstrain]';
sigmatot=sigma';
outdata=table(etot,sigmatot,p,f,dei,compS);
writetable(outdata,[filename,'.xlsx'],'Sheet',1,'Range','A1')
writetable(constsFIT,[filename,'.xlsx'],'Sheet',2,'Range','A1')

filename=strcat('Summary Data from fitted ',sampname,' data WITH NO twinning');
settings2=settings;
settings2.twipmodel=0;
settings2.plotq=0;
[sigma,p,f,dei,~,compS]=f_constituativemodel(constsFIT,exp_dat,user_dat,settings2);
p_comparison(modelstrain,compS,f,sigma,stress_dat,user_dat,constsFIT,settings2,exp_dat)
h=figure(3);
saveas(h,filename,'tiffn')
saveas(h,filename,'fig')
etot=[user_dat.modelstrain]';
sigmatot=sigma';
outdata=table(etot,sigmatot,p,f,dei,compS);
writetable(outdata,[filename,'.xlsx'],'Sheet',1,'Range','A1')
writetable(constsFIT,[filename,'.xlsx'],'Sheet',2,'Range','A1')
