%% ===================================================================== %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % Two goals of this scripts (and priority)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       1) Check the VIF values (minor)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %       2) make contrast images (minor)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %           2-1) Six categories versus base line (There are some subject
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %           data without six categoreis, i.e., four or six categories) -?
% ----------------------------------------------------------------------  %
%   * START * 
%       3) Train classifier     (major)
% ----------------------------------------------------------------------  %
%           3-1) Train data = Simple movement run No.02  (Run 06/ N = 55)
%               3-1-1) One versus others (Multi-class / Balanced design)
%               3-1-2) Test various Ridge parameter 
%               3-1-3) Linear Support Vector Machine 
%               3-1-4) y: six categories of Target angle information 
%                      x: six contrast (vs. baseline) map 
%           3-2) Test data = Simple movement run No.01  (Run 02/ N = 41)
% ====================================================================== %%

% Written by Suhwan Gim 
% (2019. 01. 22) 

% :: Examples of 'predict' function
%
% [cverr, stats, optout] = ...
%           fmri_data.predict(dat, 'algorithm_name', 'cv_svm','nfolds', 5,...
%           'MultiClass', 'error_type', 'mse','Balanced',0.01);

%
% pred_model = {};
% 
% for model_i = md
%     
%     wh_fold = [];
%     
%     dat = fmri_data;
%     dat.dat = [];
%     
%     for subj_i = 1:nsubj
%         wh_fold = [wh_fold; repmat(subj_i, 2, 1)];
%         dat.dat = [dat.dat, flat_cont_r{subj_i}{model_i}(:), flat_caps_r{subj_i}{model_i}(:)]; % concatenate times in row-wise
%         dat.Y = [dat.Y; resting.allint(subj_i); capsaicin.allint(subj_i)];
%     end
%     
%     [pred_model{model_i}.cverr, pred_model{model_i}.stats, pred_model{model_i}.optout] = fmri_data.predict(dat, 'algorithm_name', 'cv_pcr', 'nfolds', wh_fold, 'error_type', 'mse');
%     
% end
%% SETUP: directory
basedir = '/Volumes/sein/hbmnas';
scriptdir = '/Volumes/sein/hbmnas/project/SEMIC/scripts/imaging';
img_dir = fullfile(basedir, '/data/SEMIC/imaging/preprocessed');
% Model directory name
modelNameList{1} = 'model90_RUN01_Simple_motor_Single_trial';
modelNameList{2} = 'model90_RUN02_Simple_motor_Single_trial';

%  - Angle information of target dot  (5~6 Conds: 1 to 6)

subjects = canlab_list_subjects(img_dir, 'sub-semic*');


%% Train classifier 
pred_model = {};
contrasts_name{1} = '(deg_)';
contrasts_name{2} = '(Level_1)'; % (categorized)
contrasts_name{3} = '(Level_2)';
contrasts_name{4} = '(Level_3)';
contrasts_name{5} = '(Level_4)';
contrasts_name{6} = '(Level_5)';
contrasts_name{7} = '(Level_6)';

do_mask = false;
% --- NOTES ------------------------------------------------------------
% 'MultiClass' as optional input.  Important - Obj.Y must be a matrix (data x
%         class) with a column of 1 and -1 indicating each
%         class.  For example, if using 3 classes, then obj.Y
%         must have 3 columns.
% ----------------------------------------------------------------------
bal_dir = [0.01];
for bal_i = 1:2
    
    model_i=1;
    modelName = modelNameList{model_i};
    modeldir = fullfile(basedir, 'project/SEMIC/analysis/imaging/first_level',modelName); %SEE below
    outputbasedir = fullfile('/Volumes/sein/hbmnas/project/SEMIC/data/model',modelName);
    if ~isdir(outputbasedir); mkdir(outputbasedir); end
    idx = [];
    [~,a]=fileparts(modeldir);
    if strcmp(a,'model90_RUN01_Simple_motor_Single_trial') % for trainging classifier
        idx = 19:59;
    elseif strcmp(a,'model90_RUN02_Simple_motor_Single_trial') % for testing classifier
        idx = [4:20, 22:59]; % 21 skipped
        %idx = 4:59;
    end
    
    wh_fold = [];
    
    train = fmri_data;
    train.dat = [];
    train.Y = [];
    
    for cont_i = 2:length(contrasts_name) % start from 2nd contrast Level_1
        clear dat 
        load(fullfile(outputbasedir,contrasts_name{cont_i},'data_obj.mat'),'dat');
        % = apply_mask(dat,which('gray_matter_mask.nii')); 
        y_ones = ones(1,6).*(-1);        
        y_ones(cont_i-1) = 1;
        
        % make somatomotor mask
        if bal_i==2            
            gm_mask=fmri_data(which('gray_matter_mask.img'));
            Fan_Nine_Network=fmri_data(which('Fan_et_al_atlas_r280_nine_networks.nii'));
            load(which('cluster_Fan_Net_r280.mat'),'cluster_Fan_Net');
            somatomotor_mask=region2imagevec(cluster_Fan_Net.r_2mm(cluster_Fan_Net.dat(:,3) ~= 2));
            resam_somatomotor_mask=resample_space(somatomotor_mask,which('gray_matter_mask.img'),'nearest');
            gm_mask.dat(resam_somatomotor_mask.dat ~= 0) = 0;
            %apply mask
            dat = apply_mask(dat,gm_mask);
        end
        
        for image_i = 1 : length(dat.image_names)
            [~,b]=fileparts(fileparts(dat.fullpath(image_i,:)));            
            idx_num=regexp(b,'\d');
            sub_i = str2num(b(idx_num));
            wh_fold = [wh_fold, repmat(sub_i,1,1)];            
            train.dat = [train.dat, dat.dat(:,image_i)];
            train.Y = [train.Y; y_ones];
        end             
    end
    
    [pred_model{bal_i}.cverr, pred_model{bal_i}.stats, pred_model{bal_i}.optout] = ...
              predict(train, 'algorithm_name', 'cv_svm','nfolds', wh_fold,...
              'MultiClass', 'Balanced',bal_dir(1), 'error_type', 'mcr');
    
end
%% Plot 
[freq_matrix, perc_matrix] = SEMIC_confusion_matrix(pred_model{1}.stats.Y, pred_model{1}.stats.yfit,'plot');
[freq_matrix, perc_matrix] = SEMIC_confusion_matrix(pred_model{2}.stats.Y, pred_model{2}.stats.yfit,'plot');
%% Save structures
obj.Descritions = ...   
       {'';...
           '3) Train classifier     (major)'; ...
 '---------------------------------------------------------------------- ';...
 '         3-1) Train data = Simple movement run No.02  (Run 06/ N = 55) '; ...
 '              3-1-1) One versus others (Multi-class / Balanced design) ';...
'               3-1-2) Test various Ridge parameter ';...
 '              3-1-3) Linear Support Vector Machine ';...
  '             3-1-4) y: six categories of Target angle information ';...
   '                   x: six contrast (vs. baseline) map '; ...
    '       3-2) Test data = Simple movement run No.01  (Run 02/ N = 41)';...
 '====================================================================== ';...
 '';};
obj.ScriptsFullpath = '/Volumes/sein/hbmnas/project/SEMIC/scripts/3_analysis/imaging/SEMIC_190123_CV_SVM_Movement_weight_map.m';
obj.ContrastName = contrasts_name;
obj.pred_model = pred_model

save(fullfile('/Volumes/sein/hbmnas/project/SEMIC/data/model/SVM_01_simple_movement_task','SVM_01_results.mat'),'obj','-v7.3');
%%
sum(pred_model{1}.stats.Y == pred_model{1}.stats.yfit)/size((pred_model{1}.stats.yfit),1)

[freq,pro]=SEMIC_confusion_matrix(pred_model{1}.stats.yfit,pred_model{1}.stats.Y);
pred_model{1}.stats.Y 
pred_model{1}.stats.yfit

err = obj.Y ~= round(yfit); %10/7/12: Luke Chang: this will allow mcr to also be calculated using predicted probabilities

%cverr = sum(err) ./ length(err);
cverr = sum(err) ./ length(err);

phi = corr(obj.Y, yfit); %10/7/12: Luke Chang: this will calculate phi correlation coefficient between two binary variables
%orthviews(pred_model{model_i}.stats.weight_obj,which('gray_matter_mask.img'))

%% apply mask 
% applying (test data)
% (.dat .* .dat)

pred_model{1}.stats.weight_obj.volInfo;
% Weight_map 1 to 6: six catergory of angle movement
model_i = 1;
modelName = modelNameList{model_i};
modeldir = fullfile(basedir, 'project/SEMIC/analysis/imaging/first_level',modelName); %SEE below
outputbasedir = fullfile('/Volumes/sein/hbmnas/project/SEMIC/data/model',modelName);
if ~isdir(outputbasedir); mkdir(outputbasedir); end

[~,a]=fileparts(modeldir);

%%
load(fullfile(outputbasedir,contrasts_name{2},'data_obj.mat'),'dat');
ref_dat=dat;
idx = 19:59; %for run no.1
for test_i = 2:7 % index for test dat 
    clear dat weight_map
    load(fullfile(outputbasedir,contrasts_name{test_i},'data_obj.mat'),'dat');
    
    for cont_i = 2:length(contrasts_name) % start from 2nd contrast Level_1
        % weitgt_map
        cont_c = cont_i-1;        
        weigth_map = fmri_data;
        for w_i = 1:2            
            if w_i==2
                gm_mask=fmri_data(which('gray_matter_mask.img'));
                Fan_Nine_Network=fmri_data(which('Fan_et_al_atlas_r280_nine_networks.nii'));
                load(which('cluster_Fan_Net_r280.mat'),'cluster_Fan_Net');
                somatomotor_mask=region2imagevec(cluster_Fan_Net.r_2mm(cluster_Fan_Net.dat(:,3) ~= 2));
                resam_somatomotor_mask=resample_space(somatomotor_mask,which('gray_matter_mask.img'),'nearest');
                gm_mask.dat(resam_somatomotor_mask.dat ~= 0) = 0;
                %apply mask
                dat = apply_mask(dat,gm_mask);    
                
                weight_map = apply_mask(ref_dat, gm_mask);
            end            
            
            weigth_map.dat=pred_model{w_i}.stats.weight_obj.dat(:,cont_c); % 1-6
            weigth_map.volInfo = ref_dat.volInfo;
            % apply_mask
            pattern_exp_values{w_i}{test_i-1}(cont_c,:) = apply_mask(dat, weigth_map, 'pattern_expression', 'ignore_missing');
            
        end
        
        %         for i=1:6
        %             [pattern_exp_values] = apply_mask(dat, pred_model{model_i}.stats.weight_obj, 'pattern_expression', 'ignore_missing');
        %             dot_pr{cont_c}(sub_i,i) = dat.dat(:,sub_i)' * weight_map(:,i);
        %             %dot_pr{cont_c}(sub_i,i) = dat.dat(:,sub_i)' * weight_map(:,i);
        %             cor_coeff{cont_c}(sub_i,i) = corr(dat.dat(:,sub_i), weight_map(:,i),'rows','pairwise');
        %         end
        %dot_pr{cont_c}(sub_i,:)
    end
end
%%
%dot_pr{cont_c}(sub_i,:)
%cont_i = 
%dot_pr{4} == max(dot_pr{4},[],2)
for i=1:6
    for w_i=1:2
        [R{w_i},C{w_i}] =find( pattern_exp_values{w_i}{i}' == max(pattern_exp_values{w_i}{i}',[],2) );
        for ii=1:6
            idx=find(C{w_i} == ii);
            test_perc{w_i}(i,ii)=length(find(C{w_i} == ii))./size([R{w_i},C{w_i}],1) ;
        end
    end
       
end
test_perc=test_perc';
%%
c_map_str= 'gray';
subplot(1,2,1);imagesc(perc_matrix);colormap(c_map_str);colorbar;title('training sets (n=55)');
subplot(1,2,2);imagesc(test_perc);colormap(c_map_str);colorbar;title('Test sets (n=41)');