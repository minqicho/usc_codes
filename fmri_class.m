classdef fmri_class < handle
    
    properties
        number_of_subjects ;
        cortical_surface ;
        number_of_vertices ;
        number_of_faces ;
        surface_connectivity_matrix ;
        surface_connectivity_list ;
        fmri_time_series ;
        number_of_time_samples ;
        clusters_labels ;
        number_of_clusters ;
        clusters_vertices_lists ;
        clusters_average_time_series ;
        clusters_correlation_matrices ;
        clusters_inverse_correlation_matrices ;
        boundary_vertices ;
        boundary_faces ;
        boundary_rois ;
        history ;
    end;
    
    methods
        
        function self = fmri_class()
            self.number_of_subjects = 0 ;
            self.number_of_vertices = [] ;
            self.fmri_time_series = {} ;
            self.clusters_labels = {} ;
            self.clusters_vertices_lists = {};
            self.clusters_average_time_series = {};
            self.clusters_correlation_matrices = {};
            self.clusters_inverse_correlation_matrices = {};
        end;
        
        function data_validity = check_data_validity(self, number_of_vertices_local)
            data_validity = false ;
            if isempty(self.number_of_vertices)
                data_validity = true ;
                self.number_of_vertices = number_of_vertices_local ;
            elseif self.number_of_vertices == number_of_vertices_local
                data_validity = ture ;
            end;
        end;
        
        function set_cortical_surface(self, cortical_surface_local)
            [number_of_vertices_local, second_dimension] = size(cortical_surface_local.vertices) ;
            if self.check_data_validity(number_of_vertices_local)
                self.cortical_surface = cortical_surface_local ;
                self.number_of_faces = size(self.cortical_surface.faces, 1) ;
            elseif self.check_data_validity(second_dimension)
                self.cortical_surface = cortical_surface_local ;
                self.cortical_surface.vertices = self.cortical_surface.vertices' ;
                self.cortical_surface.faces = self.cortical_surface.faces' ;
                self.number_of_faces = size(self.cortical_surface.faces, 1) ;
            else
                error('Number of vertices does not match')
            end;
            [self.surface_connectivity_list, self.surface_connectivity_matrix] = vertices_connectivity_fast(self.cortical_surface) ;
        end;
        
        function set_fmri_time_series(self, fmri_time_series_local, subject_id)
            [number_of_vertices_local, number_of_time_samples_local] = size(fmri_time_series_local) ;
            if subject_id > self.number_of_subjects
                self.number_of_subjects = self.number_of_subjects + 1 ;
                subject_id = self.number_of_subjects ;
            end;
            if self.check_data_validity(number_of_vertices_local)
                self.fmri_time_series{subject_id} = fmri_time_series_local ;
                self.number_of_time_samples = number_of_time_samples_local ;
            elseif self.check_data_validity(number_of_time_samples_local)
                self.fmri_time_series{subject_id} = fmri_time_series_local' ;
                self.number_of_time_samples = number_of_vertices_local ;
            else
                error('Number of vertices does not match')
            end;
            self.pre_process_time_series(subject_id) ;
        end;
             
        function set_clusters_labels(self, clusters_labels_local, subject_id, number_of_clusters_local)
            [number_of_vertices_local, second_dimension] = size(clusters_labels_local) ;
            if subject_id > self.number_of_subjects
                self.number_of_subjects = self.number_of_subjects + 1 ;
                subject_id = self.number_of_subjects ;
            end;
            if (self.check_data_validity(number_of_vertices_local)) && (second_dimension == 1)
                self.clusters_labels{subject_id} = clusters_labels_local ;
            elseif (self.check_data_validity(second_dimension)) && (number_of_vertices_local == 1)
                self.clusters_labels{subject_id} = clusters_labels_local' ;
            else
                error('Number of vertices does not match')
            end;
            if exist('number_of_clusters_local', 'var')
                self.number_of_clusters = number_of_clusters_local ;
            else
                self.number_of_clusters = max(self.clusters_labels) ;
            end;
        end;
           
        function pre_process_time_series(self, subject_id)
            %%
            if exist('subject_id', 'var')
                subject_list = subject_id ;
            else
                subject_list = 1 : self.number_of_subjects ;
            end;
            %%
            tmp_ones_vector = ones(1, self.number_of_time_samples) ;
            for current_subject = subject_list
                current_time_series = self.fmri_time_series{current_subject} ;
                current_time_series(isnan(current_time_series)) = 0 ;
                mean_time_series_vertices = mean(current_time_series, 2) ;
                mean_time_series_matrix = mean_time_series_vertices * tmp_ones_vector ;
                current_time_series = current_time_series - mean_time_series_matrix ;
                variance_time_series_vertices = var(current_time_series, 0, 2) ;
                variance_time_series_matrix = variance_time_series_vertices * tmp_ones_vector ;
                current_time_series = current_time_series ./ variance_time_series_matrix ;
                current_time_series(isnan(current_time_series)) = 0 ;
                self.fmri_time_series{current_subject} = current_time_series ;
            end;       
        end;
        
        function compute_clusters_vertices_lists(self, subject_id)
           %%
           if exist('subject_id', 'var')
               subject_list = subject_id ;
           else
               subject_list = 1 : self.number_of_subjects ;
               self.clusters_vertices_lists = cell(self.number_of_subjects, 1) ;
           end;
           %%
           tmp_vertices_index = 1 : self.number_of_vertices ;
           for current_subject = subject_list
               self.clusters_vertices_lists{current_subject} = cell(self.number_of_clusters, 1) ;
               for current_cluster = 1 : self.number_of_clusters
                   current_cluster_vertices_index = (self.clusters_labels{current_subject} == current_cluster) ;
                   current_cluster_vertices_list = tmp_vertices_index(current_cluster_vertices_index) ;
                   self.clusters_vertices_lists{current_subject}{current_cluster} = current_cluster_vertices_list ;
               end;
           end;
        end;
        
        function compute_clusters_average_time_series(self, subject_id)
           %%
           if exist('subject_id', 'var')
               subject_list = subject_id ;
           else
               subject_list = 1 : self.number_of_subjects ;
               self.clusters_average_time_series = cell(self.number_of_subjects, 1) ;
           end;
           %%
           for current_subject = subject_list
               tmp_average_time_series = zeros(self.number_of_clusters, self.number_of_time_samples) ;
               for current_cluster = 1 : self.number_of_clusters
                   current_cluster_time_series = self.fmri_time_series{current_subject}(self.clusters_vertices_lists, :) ;
                   tmp_average_time_series(current_cluster, :) = mean(current_cluster_time_series, 1) ;
               end;
               self.clusters_average_time_series{current_subject} = tmp_average_time_series ;
           end;
        end;
        
        function compute_clusters_correlation_matrices(self, subject_id)
           %%
           if exist('subject_id', 'var')
               subject_list = subject_id ;
           else
               subject_list = 1 : self.number_of_subjects ;
               self.clusters_correlation_matrices = cell(self.number_of_subjects, 1) ;
           end;
           %%
           for current_subject = subject_list
                tmp_correlation_matrix = cov(self.compute_clusters_average_time_series{current_subject}') ;
                self.clusters_correlation_matrices{current_subject} = tmp_correlation_matrix ;      
           end;
        end;
        
        function find_boundary_vertices(self, subject_id)
            faces_labels = self.clusters_labels{subject_id}(self.cortical_surface.faces) ;
            [C, ~, IC] = unique(faces_labels, 'rows');
            number_of_unique_faces = size(C, 1) ;
            number_of_faces_labels = zeros(number_of_unique_faces, 1) ;
            for current_unique_faces = 1 : number_of_unique_faces
                [~, IA, ~] = unique(C(current_unique_faces, :)) ;
                number_of_faces_labels(current_unique_faces) = nnz(IA) ;
            end;
            boundary_faces_local = find(number_of_faces_labels > 1) ;
            boundary_vertices_index = false(self.number_of_vertices, 1) ;
            boundary_rois_local = [] ;
            for current_face = boundary_faces_local
                boundary_rois_local = union(boundary_rois_local, C(current_face, :)) ;
                boundary_vertices_index(self.cortical_surface.faces(IC == current_face)) = true ;
            end
            boundary_vertices_index(self.clusters_labels{subject_id} == 1) = false ;
            boundary_rois_local(boundary_rois_local == 1) = [] ;
            tmp_all_vertices = 1 : self.number_of_vertices ;
            self.boundary_faces{subject_id} = boundary_faces_local ;
            self.boundary_vertices{subject_id} = tmp_all_vertices(boundary_vertices_index) ;
            self.boundary_rois{subject_id} = boundary_rois_local ;
        end;
        
        function estimate_clusters_inverse_correlation_matrices_group(self, lambda, rho, max_iterations)
            %%
            if ~ exist('max_iterations', 'var')
                max_iterations = 100 ;
            end;
            self.clusters_inverse_correlation_matrices = cell(self.number_of_subjects, 1) ;
            %%
            tmp_all_correlation_matrices = zeros(self.number_of_clusters, self.number_of_clusters, self.number_of_subjects) ;
            tmp_all_inverse_correlation_matrices = zeros(size(tmp_all_correlation_matrices)) ;
            for current_subject = 1 : self.number_of_subjects
                tmp_all_correlation_matrices(:, :, current_subject) = self.clusters_correlation_matrices{current_subject} ;
                if cond(tmp_all_correlation_matrices(:, :, current_subject)) < 1e5
                    tmp_all_inverse_correlation_matrices(:, :, current_subject) = inv(tmp_all_correlation_matrices(:, :, current_subject)) ;
                else
                    tmp_all_inverse_correlation_matrices(:, :, current_subject) = pinv(tmp_all_correlation_matrices(:, :, current_subject)) ;
                end;
            end;
            tmp_all_inverse_correlation_matrices(isnan(tmp_all_inverse_correlation_matrices)) = 0 ;
            g_matrix = zeros(size(tmp_all_inverse_correlation_matrices)) ;
            u_matrix = g_matrix - tmp_all_inverse_correlation_matrices ;
            theta = lambda / rho ;
            for current_iteration = 1 : max_iterations
                for current_subject = 1 : self.number_of_subjects
                [q_matrix, l_matrix] = eig(rho * (g_matrix(:, :, current_subject) - u_matrix(:, :, current_subject))-tmp_all_correlation_matrices(:, :, current_subject)) ;
                a_matrix = diag((diag(l_matrix) + sqrt(diag(l_matrix) .^ 2 + 4 * rho)) ./ 2 ./ rho) ;
                tmp_all_inverse_correlation_matrices(:, :, current_subject) = q_matrix * a_matrix * q_matrix ;
                end;
                g_vector = zeros(self.number_of_clusters * self.number_of_clusters, self.number_of_subjects) ;
                u_vector = reshape(u_matrix, self.number_of_clusters * self.number_of_clusters, self.number_of_subjects) ;
                tmp_all_inverse_correlation_vector = reshape(tmp_all_inverse_correlation_matrices, self.number_of_clusters * self.number_of_clusters, self.number_of_subjects) ;
                for current_clusters_pair = 1 : self.number_of_clusters * self.number_of_clusters
                    current_clusters_pair_norm = norm(tmp_all_inverse_correlation_vector(current_clusters_pair, :) + u_vector(current_clusters_pair, :)) ;
                    if current_clusters_pair_norm > theta
                        g_vector(current_clusters_pair, :) = (current_clusters_pair_norm - theta) / current_clusters_pair_norm * (tmp_all_inverse_correlation_vector(current_clusters_pair, :) + u_vector(current_clusters_pair, :)) ;
                    end;
                end;
                g_matrix = reshape(g_vector, self.number_of_clusters, self.number_of_clusters, self.number_of_subjects) ;
                u_matrix = u_matrix + rho * (tmp_all_inverse_correlation_matrices - g_matrix) ;
            end;
            tmp_all_inverse_correlation_matrices = g_matrix ;
            for current_subject = 1 : self.number_of_subjects
                self.clusters_inverse_correlation_matrices{current_subject} = tmp_all_inverse_correlation_matrices(:, :, current_subject) ;
            end;
        end;
        
        function similarity_matrix = compute_similarity_matrix(self, subject_id, lambda)
            number_of_boundary_rois = length(self.boundary_rois{subject_id}) ;
            number_of_boundary_vertices = length(self.boundary_vertices{subject_id});
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
                current_vertex = self.boundary_vertices{subject_id}(current_vertex_index);
                vertex_time_series = self.fmri_time_series{subject_id}(current_vertex, :) ;
                neighbor_rois = union(self.clusters_labels{subject_id}(current_vertex), self.clusters_labels{subject_id}(self.surface_connectivity_list(current_vertex))) ;
                neighbor_rois = intersect(neighbor_rois, self.boundary_rois{subject_id}) ;
                neighbor_rois_index = zeros(size(neighbor_rois)) ;
                for current_roi_index = 1 : length(neighbor_rois)
                    neighbor_rois_index(current_roi_index) = find(self.boundary_rois{subject_id} == neighbor_rois(current_roi_index));
                end;
                neighbor_rois_index(neighbor_rois_index == 0) = [] ;
                for current_roi_index = neighbor_rois_index
                    current_roi = self.boundary_rois{subject_id}(current_roi_index) ;
                    roi_all_time_series = self.clusters_average_time_series{subject_id} ;
                    roi_inverse_correlation_vector = self.clusters_inverse_correlation_matrices{subject_id}(:, current_roi) ;
                    roi_vertices_count = nnz(self.clusters_vertices_lists{subject_id}{current_roi}) ;
                    tmp_similarity = (2 * vertex_time_series * roi_all_time_series' * roi_inverse_correlation_vector) / roi_vertices_count;
                    if self.clusters_labels{subject_id}(current_vertex) == current_roi
                        tmp_similarity = tmp_similarity + 2 * lambda ;
                    end;
                    tmp_similarity = exp(tmp_similarity) ;
                    similarity_matrix(current_roi_index, current_vertex_index + number_of_boundary_rois) = tmp_similarity ;
                    similarity_matrix(current_vertex_index + number_of_boundary_rois, current_roi_index) = similarity_matrix(current_roi_index, current_vertex_index + number_of_boundary_rois) ;
                end;
                for neighbor_vertex = self.surface_connectivity_list(current_vertex)
                    [is_neighbor, neighbor_vertex_index] = ismember(neighbor_vertex, self.boundary_vertices{subject_id}) ;
                    if is_neighbor
                        similarity_matrix(current_vertex_index + number_of_boundary_rois, neighbor_vertex_index + number_of_boundary_rois) = exp(lambda) ;
                        similarity_matrix(neighbor_vertex_index + number_of_boundary_rois, current_vertex_index + number_of_boundary_rois) = similarity_matrix(current_vertex_index + number_of_boundary_rois, neighbor_vertex_index + number_of_boundary_rois) ;
                    end;
                end;
            end;
            similarity_matrix(similarity_matrix == 1) = 0 ;
        end;
        
        function relabel_boundary_vertices_from_graph_cut_results(self, subject_id, graph_cut_results)
            number_of_boundary_rois = length(self.boundary_rois{subject_id}) ;
            number_of_boundary_vertices = length(self.boundary_vertices{subject_id});
            rois_results = graph_cut_results(1 : number_of_boundary_rois, :) ;
            vertices_results = graph_cut_results(number_of_boundary_rois + 1 : number_of_boundary_rois + number_of_boundary_vertices, :) ;
            rois_clusters = rois_results * (1 : size(rois_results, 2)) ;
            vertices_clusters = vertices_results * (1 : size(vertices_results, 2)) ;
            for current_roi = 1 : number_of_boundary_rois
               candidate_vertices = self.boundary_vertices(vertices_clusters == rois_clusters(current_roi)) ;
               self.clusters_labels{subject_id}(candidate_vertices) = self.boundary_rois(current_roi) ;
            end
        end
        
        function estimate_parcellation(self, subject_id, lambda)
            self.find_boundary_vertices(subject_id) ;
            similarity_matrix = self.compute_similarity_matrix(subject_id, lambda) ;
            graph_cut_results = ncutW(similarity_matrix, length(self.boundary_rois{subject_id})) ;
            self.relabel_boundary_vertices_from_graph_cut_results(subject_id, graph_cut_results) ;
        end;
              
        function minqi_group_clustering(self, max_iterations, lambda1, lambda2)
            %%
            if isempty(self.cortical_surface)
                error('Need to have a cortical surface!') ;
            end;
            if isempty(self.fmri_time_series)
                error('Need to have fMRI time series!') ;
            end;
            if isempty(self.clusters_labels)
                error('Need to have an initialize parcellation!') ;
            end;
            if ~ exist('max_iterations', 'var')
                max_iterations = 20 ;
            end;
            if ~ exist('lambda1', 'var')
                lambda1 = 1 ;
            end;
            if ~ exist('lambda2', 'var')
                lambda2 = 0.2 ;
            end;
            %%
            for current_iteration = 1 : max_iterations
                self.compute_clusters_vertices_lists ;
                self.compute_clusters_average_time_series ; 
                self.compute_clusters_correlation_matrices ;
                self.estimate_clusters_inverse_correlation_matrices_group(lambda2, 1) ;
                tmp_parcellations = self.clusters_labels ;
                self.boundary_vertices = cell(self.number_of_subjects, 1);
                self.boundary_faces = cell(self.number_of_subjects, 1);
                self.boundary_rois = cell(self.number_of_subjects, 1);
                parfor current_subject = 1 : self.number_of_subjects
                    tmp_parcellations{current_subject} = self.estimate_parcellation(current_subject, lambda1) ;
                end;
                self.clusters_labels = tmp_parcellations ;
            end;
            
        end;
            
    end;
end
