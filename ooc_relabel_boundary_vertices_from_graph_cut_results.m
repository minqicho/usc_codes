function refined_parcellation = ooc_relabel_boundary_vertices_from_graph_cut_results( subject_instance, graph_cut_results, subject_similarity_matrix)
%OOC_RELABEL_BOUNDARY_VERTICES_FROM_GRAPH_CUT_RESULTS Summary of this function goes here
%   Detailed explanation goes here
number_of_boundary_rois = length(subject_instance.boundary_rois) ;
number_of_boundary_vertices = length(subject_instance.boundary_vertices);
rois_results = graph_cut_results(1 : number_of_boundary_rois, :) ;
vertices_results = graph_cut_results(number_of_boundary_rois + 1 : number_of_boundary_rois + number_of_boundary_vertices, :) ;
rois_clusters = rois_results * (1 : size(rois_results, 2))' ;
vertices_clusters = vertices_results * (1 : size(vertices_results, 2))' ;
for current_roi_index = 1 : number_of_boundary_rois
    current_roi = subject_instance.boundary_rois(current_roi_index) ;
    candidate_vertices = subject_instance.boundary_vertices(vertices_clusters == rois_clusters(current_roi_index)) ;
    subject_instance.clusters_labels(candidate_vertices) = current_roi ;
end
%%
rois_clusters_sorted = sort(rois_clusters) ;
duplicate_clusters = find(diff(rois_clusters_sorted) == 0) ;
if ~isempty(duplicate_clusters)
    for current_duplicate_cluster_index = 1 : length(duplicate_clusters)
        current_duplicate_clusters = find(rois_clusters == rois_clusters_sorted(duplicate_clusters(current_duplicate_cluster_index))) ;
        current_duplicate_vertices = (vertices_clusters == rois_clusters_sorted(duplicate_clusters(current_duplicate_cluster_index))) ;
        duplicate_similarity_matrix = subject_similarity_matrix(current_duplicate_clusters, [false(number_of_boundary_rois, 1) ; current_duplicate_vertices]);
        [~, max_roi_indices] = max(duplicate_similarity_matrix, [], 1);
        subject_instance.clusters_labels(subject_instance.boundary_vertices(current_duplicate_vertices)) = subject_instance.boundary_rois(current_duplicate_clusters(max_roi_indices)) ;       
    end
end
%%
refined_parcellation = subject_instance.clusters_labels ;
end

