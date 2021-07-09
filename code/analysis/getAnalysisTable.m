
function analysis_table = getAnalysisTable(study_settings)
    analysis_table_body = study_settings.get('analysis_table', 1);
    analysis_table_header = study_settings.get('analysis_table_header', 1);
    
    analysis_table = cell2table(analysis_table_body);
    analysis_table.Properties.VariableNames = analysis_table_header;
end
