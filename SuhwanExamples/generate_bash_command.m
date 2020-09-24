
jobsn = 3;
%% 
[dir,str] =  SEMIC_HPC_set_dir('HPC');
%%
home_dir = '/home/suhwan/Desktop';
output_home_dir = 'Matlab_output/';
log_folder = fullfile(home_dir, output_home_dir);
%%
clc
str = [];
txtfile=fopen(fullfile(home_dir,[date '_generate_bash.txt']),'w');
for tr_i = 1:45
    %outputdir = fullfile(dir.prj_dir,'/scripts/4_custom_function/Mediation_dream_test','output',sprintf('seg_%02d',seg_i));
    outputdir = fullfile('/sas1/cocoanlab/data/SEMIC/190309_Mediation_dream_test/output',sprintf('tr_%02d',tr_i));
    for sec_i = 1:32
    for job_i = 1:jobsn 
        
        code_filename = fullfile(outputdir, [sprintf('SEMIC_mediation_brain_%02d_of_32_run_suhwan',sec_i) '_' sprintf('%03d',job_i) '.m']);
        
        scripts_path = code_filename;
        str = ['nohup matlab_orig -nodesktop -nodisplay -nosplash -singleCompThread -r " run ' scripts_path ' " > output ' fullfile(log_folder,sprintf('_%03d.txt',job_i)) ' < /dev/null';];
        str2 = [str '\r\n'];
        fprintf(txtfile,str2);
        disp(str);
    end
    end
end
%%
clc
jobsn = 3;
str = [];
txtfile=fopen(fullfile(home_dir,[date '_generate_bash.txt']),'w');
for tr_i = 1:45
    %outputdir = fullfile(dir.prj_dir,'/scripts/4_custom_function/Mediation_dream_test','output',sprintf('seg_%02d',seg_i));
    outputdir = fullfile('/sas1/cocoanlab/data/SEMIC/190420_4mm_whold_brain_mediation_dream/output',sprintf('tr_%02d',tr_i));
    sec_i=32;
    for job_i = 1:jobsn
        
        code_filename = fullfile(outputdir, [sprintf('SEMIC_mediation_brain_%02d_of_32_run_suhwan',sec_i) '_' sprintf('%03d',job_i) '.m']);
        
        scripts_path = code_filename;
        str = ['nohup matlab_orig -nodesktop -nodisplay -nosplash -singleCompThread -r " run (''' scripts_path '''); quit " > output ' fullfile(log_folder,sprintf('_%03d.txt',job_i)) ' < /dev/null';];
        str2 = [str '\r\n'];
        fprintf(txtfile,str2);
        disp(str);
    end
    
end
%%
%for i=1:numel(str), eval(str{i}); end;
str = ['no hub matlab_orig -nodesktop -nodisplay -nosplash -r " ' scripts_path ' " > output.txt < /dev/null';];


% run_preproc_command = ['matlab -nodesktop -nosplash -nodisplay -r "addpath(genpath(''' fullfile(rscdir, 'matlab_toolboxes/surface_preprocessing') '''));' ...
%     'CAPS2_preproc([' num2str(sj_num) ']); quit" >> ~/' sprintf('sub-caps%.3d_log.txt', sj_num) ' 2>&1 < /dev/null &'];
clipboard('copy', str_preproc_command);

%% at once?
sj_num = 1:53;
run_preproc_command = [];
for i = 1:numel(sj_num)
run_preproc_command = [run_preproc_command 'matlab_orig -nodesktop -nosplash -nodisplay -r '  ...
 '"run(''' fullfile(basedir, ['projects/CAPS2/scripts/onlyfan([' num2str(sj_num(i)) '])']) '''); quit"'  ...
 ' >> ~/' sprintf('sub-caps%.3d_temp_log.txt', sj_num(i)) ' 2>&1 < /dev/null &'];
end
clipboard('copy', run_preproc_command);

%%
run_preproc_command = [];
for i = 1:45
    for job_i = 1:jobsn
        
        code_filename = fullfile(outputdir, [sprintf('SEMIC_mediation_brain_%02d_of_32_run_suhwan',sec_i) '_' sprintf('%03d',job_i) '.m']);
        scripts_path = code_filename;
        
        run_preproc_command = [run_preproc_command 'matlab_orig -nodesktop -nosplash -nodisplay -r '  ...
            '"run(''' scripts_path '''); quit"'  ...
            ' >> ~/' fullfile(log_folder,sprintf('_%03d.txt',job_i)) ' 2>&1 < /dev/null &'];
    end
end

%%
% process kill [,pid] 

cmdd = [];

cmdd = [cmdd 'nohup matlab_orig -nodesktop -nodisplay -nosplash -r " ' scripts_path ' " 1> output.txt  2> errorlog.txt &';];