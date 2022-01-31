function decoding_error_controls(current_folder,method)
% options are:
% 'global_remap'
% 'rate_remp'

if exist('lap_times.mat')~=2
    extract_laps;
end

switch method
    case 'global_remap'
        copyfile('lap_times.mat',['..\CONTROLS\global_remapped\' current_folder]);
        cd(['..\CONTROLS\global_remapped\' current_folder]); 
        bayesian_decoding_error('method','leave one out','control','global_remap');
        
        
    case 'rate_remap'
        copyfile('lap_times.mat',['..\CONTROLS\rate_remapped\' current_folder]);
        cd(['..\CONTROLS\rate_remapped\' current_folder]); 
        bayesian_decoding_error('method','leave one out','control','rate_remap');
end
        

end