function refined_parcellation = ooc_estimate_parcellation(subject_instance, lambda)
%OOC_ESTIMATE_PARCELLATION Summary of this function goes here
%   Detailed explanation goes here
subject_similarity_matrix = ooc_compute_similarity_matrix(subject_instance, lambda) ;
subject_graph_cut_result = ncutW(subject_similarity_matrix, subject_instance.number_of_clusters) ;
refined_parcellation = ooc_relabel_boundary_vertices_from_graph_cut_results(subject_instance, subject_graph_cut_result, subject_similarity_matrix) ;
end

