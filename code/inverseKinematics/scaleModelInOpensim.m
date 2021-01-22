%     This file is part of the CoBaL code base
%     Copyright (C) 2017 Hendrik Reimann <hendrikreimann@gmail.com>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function scaleModelInOpensim()
    % load settings
    subject_settings = SettingsCustodian('subjectSettings.txt');
    subject_settings = loadSettingsFromFile('subject');
    collection_date = subject_settings.get('collection_date');
    subject_id = subject_settings.get('subject_id')
    
    %% set up
    import org.opensim.modeling.*
    generic_model_file = 'CoBaLWalker50.osim';
    output_model_file = makeFileName(collection_date, subject_id, '.osim');
    % model is expected relative to the current directory, so make a temporary copy
    copyfile([getCobalPath filesep 'resources' filesep 'opensim' filesep generic_model_file], ['opensim' filesep generic_model_file])
    
    % setup model scaling 
    static_trial_type = subject_settings.get('static_reference_trial_type');
    static_trial_number = subject_settings.get('static_reference_trial_number');
    static_file_name = ['marker' filesep makeFileName(collection_date, subject_id, static_trial_type, static_trial_number, 'marker.trc')];
    generic_setup_file_scale = [getCobalPath filesep 'resources' filesep 'opensim' filesep 'CoBaLWalker50_setupScale.xml'];
    DOMnode = xmlread(generic_setup_file_scale);
    document_node = DOMnode.getFirstChild;
    
    % set name of resulting model
    scale_tool_node = document_node.getElementsByTagName("ScaleTool").item(0);
    scale_tool_node.setAttribute("name", output_model_file)
    
    % set generic model
    model_maker_node = scale_tool_node.getElementsByTagName("GenericModelMaker").item(0);
    model_file_node = model_maker_node.getElementsByTagName("model_file").item(0);
    child = model_file_node.getFirstChild;
    child.setData(generic_model_file)

    % set marker file name for model scaler
    model_scaler_node = scale_tool_node.getElementsByTagName("ModelScaler").item(0);
    marker_file_node = model_scaler_node.getElementsByTagName("marker_file").item(0);
    child = marker_file_node.getFirstChild;
    child.setData(static_file_name)

    % set marker file name for marker placer
    marker_placer_node = scale_tool_node.getElementsByTagName("MarkerPlacer").item(0);
    marker_file_node = marker_placer_node.getElementsByTagName("marker_file").item(0);
    child = marker_file_node.getFirstChild;
    child.setData(static_file_name)

    % set output file name for marker placer
    marker_placer_node = scale_tool_node.getElementsByTagName("MarkerPlacer").item(0);
    output_model_file_node = marker_placer_node.getElementsByTagName("output_model_file").item(0);
    child = output_model_file_node.getFirstChild;
    child.setData(output_model_file)

    % save subject-specific scale file
    setup_file_scale = [pwd filesep 'opensim' filesep makeFileName(collection_date, subject_id, 'setupScale.xml')];
    xmlwrite(setup_file_scale, DOMnode);    
    
    % scale
    scale_tool = ScaleTool(setup_file_scale);
    scale_tool.setPrintResultFiles(1)
    scale_tool.run();
    
    % remove local copy of the generic model
    delete(['opensim' filesep generic_model_file]);
    
end













