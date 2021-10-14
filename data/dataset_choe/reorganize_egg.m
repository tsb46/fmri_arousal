sess_n = 19
run = 1

run_list_n = dlmread("run_list_choe.csv", ",",1,0);
egg_struct = load('egg/raw/EGGGICA_prepro_ds_all.mat');

for i = 1:length(run_list_n)
	egg_sess_run = egg_struct.EGGGICA_prepro_ds_all.sess(i).run(1).EGG;
	dlmwrite(['egg/raw/0' num2str(run_list_n(i)) '_run1_EGG.txt'], egg_sess_run); 
end