clc
clear
%%               MMPDNB
%                  Multi-modal optimization to identify personalized biomarkers for detecting early warning signal of individual patients in cancer
%
%                  Attention : Please raad 'ReadMe' file before running MMPDNB_code     
%
%                  input  :
%                         cancer_data     :    BRCA or LUNG
%                         path            :    The path of the user where the 'Main,m' is located
%
%                   output :         
%                          PGIN_'cancer_data'      : patient of 'cancer_data' samples' personalized gene interaction network
%                          'cancer_data'_result    : non dominated solution of patient samples obtained by MMPDNB
%                          'cancer_data'_DNB_name  : gene name of MMPDNB or PDNB of 'cancer_data' patients
%
%                   Proposed by J.J.Liang, Wei feng, Guo, and Zongwei Li in Zhengzhou University on November 15, 2021.
%                   If any problem,pleasse contact zero999396@outlook.com for help.
                  

%% ***************************Input*************************************

path='...\';        % The path of 'Main,m' on the user's computer and '\' need reserve.
unzip('BRCA_tumor.zip');              % Take BRCA as an example
unzip('BRCA_normal.zip');
 
expression_tumor_fileName = strcat('BRCA_tumor.txt');
expression_normal_fileName = strcat('BRCA_normal.txt');

%% *******************The function of MMPDNB****************************

%% default parameters
popnum=300;
Max_CalNum=30000;
Experiment_num=30;

gen_name=mmpdnb(expression_tumor_fileName,expression_normal_fileName,path,popnum,Max_CalNum,Experiment_num);

%% *********Ouutput**********

data2=load([path,expression_tumor_fileName(1:4),'_clinical_stage_information.mat']);
output_file_name=strcat(expression_tumor_fileName(1:4),'_DNB_name.txt');
fid = fopen(output_file_name,'w');
        for j = 1:length(gen_name)
            k=1;
            fprintf(fid ,num2str(j));
            fprintf(fid ,'%s', ':');
            patient_name=cell2mat(data2.Final_Sample_name(j,1));
            fprintf(fid ,patient_name);
            fprintf(fid ,'%s', ' ');
            while k<=length(gen_name{j,1})
                if length(gen_name{j,1})>1
                    MM=strcat('multi-modal PDNB',num2str(k));
                else
                     MM=strcat('PDNB');
                end
                fprintf(fid ,MM); 
                fprintf(fid ,'%s', ': ');
                l=1;
                while l<=length(gen_name{j,1}{k,1})
                n=cell2mat(gen_name{j,1}{k,1}(l,1));
                fprintf(fid ,n);
                fprintf(fid ,'%s', '\');
                l=l+1;
                end
            fprintf(fid ,'%s', ';');
            k=k+1;
            end
            fprintf(fid ,'\n');
        end
    fclose(fid);