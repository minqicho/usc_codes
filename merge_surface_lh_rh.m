function both_hemi_data = merge_surface_lh_rh( left_hemi_data, right_hemi_data, fields_to_merge )
%COMBINE_LH_RH Summary of this function goes here
%   Detailed explanation goes here
number_of_lh_vertices = size(left_hemi_data.vertices, 1) ;
number_of_rh_vertices = size(right_hemi_data.vertices, 1) ;
both_hemi_data.vertices = [left_hemi_data.vertices ; right_hemi_data.vertices] ;
both_hemi_data.faces = [left_hemi_data.faces ; right_hemi_data.faces + number_of_lh_vertices] ;
both_hemi_data.number_of_lh_vertices = number_of_lh_vertices ;
if exist('fields_to_merge', 'var')
    for current_field_index = 1 : length(fields_to_merge)
       current_field_name = char(fields_to_merge{current_field_index}) ;
       if any(ismember((fieldnames(left_hemi_data)), current_field_name))
           if any(ismember((fieldnames(right_hemi_data)), current_field_name))
               if number_of_lh_vertices == size(left_hemi_data.(current_field_name), 1)
                   if number_of_rh_vertices == size(right_hemi_data.(current_field_name), 1)
                       both_hemi_data.(current_field_name) = [left_hemi_data.(current_field_name) ; right_hemi_data.(current_field_name)] ;
                   end
               end
           end
       end
    end
end
end

