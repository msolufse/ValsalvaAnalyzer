function [T] = DataSummary(list)
% function DataSummary
% Input: list of quantities to save to data file
% Output: excel spreadsheet with computed quantities
% Description: mat2excel Writes markers from mat files to an excel spreadsheet

cd .. 
Info = strings(length(list),5);
ModelParameters = zeros(length(list),10);

for i = 1:length(list)
    patient   = list{i};
    patient   = patient(1:end-4);
    load(strcat('WS/',patient,'_WS.mat'),'Weight','Height','Sex','Age')
    Info(i,1)    = patient;
    Info(i,2)    = Age;
    Info(i,3)    = Sex;
    Info(i,4)    = Height;
    Info(i,5)    = Weight;
    
    T_ClinicalRatios(i,:) = readtable(strcat('Markers/',patient,'_markers.xlsx'));
    
    load(strcat('Optimized/',patient,'_opt.mat'),'optpars')
    optpars = exp(optpars);
    ModelParameters(i,1)  = optpars(11);  %tau_pb estimated
    ModelParameters(i,2)  = optpars(12);  %tau_pr estimated
    ModelParameters(i,3)  = optpars(19);  %xi_w estimated
    ModelParameters(i,4)  = optpars(20);  %xi_pb patient specific
    ModelParameters(i,5)  = optpars(21);  %xi_pr patient specific
    ModelParameters(i,6)  = optpars(22);  %xi_s patient specific
    ModelParameters(i,7)  = optpars(23);  %H_I patient specific
    ModelParameters(i,8)  = optpars(24);  %H_pb estimated
    ModelParameters(i,9)  = optpars(25);  %H_pr estimated
    ModelParameters(i,10) = optpars(26);  %H_s estimated
    ModelParameters(i,11) = optpars(27);  %H_Ia estimated
    ModelParameters(i,12) = optpars(28);  %tchar patient specific
end
T_ClinicalRatios = removevars(T_ClinicalRatios,{'Var1'});
T_Info = array2table(Info, 'VariableNames', {'Patient', 'Age', 'Sex','Height','Weight'});
T_ModelParameters = array2table(ModelParameters, 'VariableNames', ...
    {'tau_pb', 'tau_pr', 'xi_w', 'xi_pb','xi_pr','xi_s','H_I','H_bp','H_pr','H_s','HIa','tchar'});
T = [T_Info, T_ClinicalRatios, T_ModelParameters];

size(T)

prompt = {'Enter data summary file name (ex: FileName.xlsx)'};
dlgtitle = 'Save Data';
dims = [1 36];
definput = {'FileName.xlsx','hsv'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

if isempty(answer) == 1
    return
end

if isfile(answer{1})
    disp('file exist')
    quest = {'Add or overwrite existing data?'};
    dlgtitle = 'Save Data';
    answer2 = questdlg(quest,dlgtitle,'Add','Overwrite','Add');
    if strcmp(answer2, 'Add') == 1
%       writetable(T,answer{1},'UseExcel', true, 'WriteMode','Append','AutoFitWidth',false)
        writetable(T,answer{1},'WriteMode','Append','AutoFitWidth',false)
    else
        writetable(T,answer{1}, 'WriteMode','Overwrite')
    end
else
    disp('file does not exist')
    writetable(T,answer{1}, 'WriteMode','Overwrite')
end

cd Core
end
