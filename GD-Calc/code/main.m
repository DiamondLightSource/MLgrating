function main
% Run one or more demo scripts (gdc_demo*.m); capture text and graphic
% (figure) output and save to results directory.
%
% Version 26-Jul-2022

exclude_files = {};
if 1
    disp('NOTE: GD-Calc version v4 requires MATLAB version R2022a or ')
    disp('later, which is not currently supported by Code Ocean. Three ')
    disp('files, gdc.m, gdc_engine.p, and gdc_eff.m, have been ')
    disp('replaced by their v3 versions to run in Code Ocean.')
    disp(' ')
    disp('gdc_demo16.m will not run under v3 because the last gdc ')
    disp('argument (all_inc_order) must be logical (true|false) in v3.')
    disp(' ')
    disp('To run version v4, do the following:')
    disp('  - Install MATLAB R2022a or later.')
    disp('  - Repace gdc.m, gdc_engine.p, and gdc_eff.m with the v4 ')
    disp('    versions (gdc_v4.m, gdc_engine_v4.p, and gdc_eff_v4.m).')
    disp('  - Remove the ''1f 1 ...'' code block in main.m.')
    disp(' ')
    exclude_files = {'gdc_demo16.m'};
end

% Select all demo scripts.
files = dir('gdc_demo*.m');
files = {files.name};
if ~isempty(exclude_files)
    files = setdiff(files,exclude_files);
end
for file = files
    file = file{1}; %#ok<FXSET>
    disp(['Running ' file ' ...'])
    file = erase(file,'.m'); %#ok<FXSET>
    fid = fopen(['../results/' file '.txt'],'w');
    str = evalc(file);
    fprintf(fid,str);
    fclose(fid);
    figs =  findobj('type','figure');
    for j = length(figs):-1:1
        figure(figs(j));
        N = get(gcf,'Number');
        fig_tif = ['../results/' file '_Fig' num2str(N) '.tif'];
        print(gcf,fig_tif,'-dtiff','-r300');
        close(gcf)
    end
end
end
