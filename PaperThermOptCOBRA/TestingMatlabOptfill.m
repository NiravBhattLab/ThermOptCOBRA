clear
clc
%% for toy model 1
model = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM1_TDb1\TM1.mat');
database = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM1_TDb1\TDb1.mat');

tic
[TICs_o1,Direction_o1] = OptFillTFP(model.model,database.model,1);
toc_o1=toc

tic
m_db = mergeTwoModels(model.model,database.model);
[TICs_t1,Direction_t1,~,~,~,nMILPS_1,nDVs_1] = ThermOptEnumMILP_count_milp_dvs(m_db);
toc_t1=toc

%% for toy model 2
model = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM2_TDb2\TM2.mat');
database = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM2_TDb2\TDb2.mat');
tic
[TICs_o2,Direction_o2] = OptFillTFP(model.model,database.model,1);
toc_o2=toc

tic
m_db = mergeTwoModels(model.model,database.model);
[TICs_t2,Direction_t2,~,~,~,nMILPS_2,nDVs_2] = ThermOptEnumMILP_count_milp_dvs(m_db);
toc_t2=toc

%% for toy model 3
model = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM3_TDb3\TM3.mat');
database = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\TM3_TDb3\TDb3.mat');
tic
[TICs_o3,Direction_o3] = OptFillTFP(model.model,database.model,1);
toc_o3=toc

tic
m_db = mergeTwoModels(model.model,database.model);
[TICs_t3,Direction_t3,~,~,~,nMILPS_3,nDVs_3] = ThermOptEnumMILP_count_milp_dvs(m_db);
toc_t3=toc

%% for ecoli model
model = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\iJR904_iAF1260\iJR904.mat');
database = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\TIC\deleteit_2\optfill_in_matlab\iJR904_iAF1260\iAF1260.mat'); 
tic
[TICs_o4,Direction_o4] = OptFillTFP(model.model,database.model,1);
toc_o4=toc

tic
m_db = mergeTwoModels(model.model,database.model);
[TICs_t4,Direction_t4,~,~,~,nMILPS_4,nDVs_4] = ThermOptEnumMILP_count_milp_dvs(m_db);
toc_t4=toc