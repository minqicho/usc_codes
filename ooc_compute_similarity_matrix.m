function similarity_matrix = ooc_compute_similarity_matrix(subject_instance, lambda)
%OOC_COMPUTE_SIMILARITY_MATRIX Summary of this function goes here
%   Detailed explanation goes here
number_of_boundary_rois = length(subject_instance.boundary_rois) ;
number_of_boundary_vertices = length(subject_instance.boundary_vertices);
number_of_nodes = number_of_boundary_rois + number_of_boundary_vertices ;
similarity_matrix = zeros(number_of_nodes, number_of_nodes);
for current_roi_1 = 1 : number_of_boundary_rois
    for current_roi_2 = 1 : number_of_boundary_rois
        if current_roi_1 ~= current_roi_2
            similarity_matrix(current_roi_1, current_roi_2) = -1 ;
            similarity_matrix(current_roi_2, current_roi_1) = similarity_matrix(current_roi_1, current_roi_2) ;
        end
    end;
end;
for current_vertex_index = 1 : number_of_boundary_vertices
    current_vertex = subject_instance.boundary_vertices(current_vertex_index);
    vertex_time_series = subject_instance.fmri_time_series(current_vertex, :) ;
    neighbor_rois = union(subject_instance.clusters_labels(current_vertex), subject_instance.clusters_labels(subject_instance.surface_connectivity_list{current_vertex})) ;
    neighbor_rois = intersect(neighbor_rois, subject_instance.boundary_rois) ;
    neighbor_rois(neighbor_rois == 0) = [];
    for current_roi_index = 1 : length(subject_instance.boundary_rois)
        current_roi = subject_instance.boundary_rois(current_roi_index) ;
        if ~isempty(neighbor_rois)
            if any(ismember(neighbor_rois, current_roi)) 
                roi_all_time_series = subject_instance.clusters_average_time_series;
                roi_inverse_correlation_vector = subject_instance.clusters_inverse_correlation_matrices(:, current_roi) ;
                roi_vertices_count = nnz(subject_instance.clusters_vertices_lists{current_roi}) ;
                tmp_similarity = trace( 2 * vertex_time_series * roi_all_time_series' * roi_inverse_correlation_vector) / roi_vertices_count;
                if subject_instance.clusters_labels(current_vertex) == current_roi
                    tmp_similarity = tmp_similarity + 2 * lambda ;
                end;
                similarity_matrix(current_roi_index, current_vertex_index + number_of_boundary_rois) = exp(tmp_similarity) ;
                similarity_matrix(current_vertex_index + number_of_boundary_rois, current_roi_index) = similarity_matrix(current_roi_index, current_vertex_index + number_of_boundary_rois) ;
           end
        end
    end;
    for neighbor_vertex_index = 1 : length(subject_instance.surface_connectivity_list{current_vertex})
        neighbor_vertex = subject_instance.surface_connectivity_list{current_vertex}(neighbor_vertex_index) ;
        is_neighbor_vertex_index = find(subject_instance.boundary_vertices == neighbor_vertex, 1) ;
        if ~ isempty(is_neighbor_vertex_index)
            similarity_matrix(current_vertex_index + number_of_boundary_rois, is_neighbor_vertex_index + number_of_boundary_rois) = exp(lambda) ;
            similarity_matrix(is_neighbor_vertex_index + number_of_boundary_rois, current_vertex_index + number_of_boundary_rois) = similarity_matrix(current_vertex_index + number_of_boundary_rois, is_neighbor_vertex_index + number_of_boundary_rois) ;
        end;
    end;
end;
%similarity_matrix = exp(similarity_matrix) ;
%similarity_matrix(similarity_matrix == 1) = 0 ;
similarity_matrix(isnan(similarity_matrix)) = 0 ;
similarity_matrix(isinf(similarity_matrix)) = 0 ;
end

