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
        reference_data ;
        task_z_scores ;
        task_names ;
        history ;
    end;
    
    methods
        
        function self = fmri_class()
            self.number_of_subjects = 0 ;
            self.number_of_vertices = 0 ;
            self.number_of_clusters = 0 ;
            self.fmri_time_series = {} ;
            self.clusters_labels = {} ;
            self.clusters_vertices_lists = {};
            self.clusters_average_time_series = {};
            self.clusters_correlation_matrices = {};
            self.clusters_inverse_correlation_matrices = {};
            self.history = {};
            self.reference_data = {} ;
            self.task_z_scores = {} ;
        end;
        
        function data_validity = check_data_validity(self, number_of_vertices_local)
            data_validity = false ;
            if (isempty(self.number_of_vertices)) || (self.number_of_vertices == 0)
                data_validity = true ;
                self.number_of_vertices = number_of_vertices_local ;
            elseif self.number_of_vertices == number_of_vertices_local
                data_validity = true ;
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
        
        function set_reference_data(self, reference_data_local, subject_id)
            [number_of_vertices_local, second_dimension] = size(reference_data_local) ;
            if self.check_data_validity(number_of_vertices_local)
                self.reference_data{subject_id} = reference_data_local ;
            elseif self.check_data_validity(second_dimension)
                self.reference_data{subject_id} = reference_data_local' ;
            else
                error('Number of vertices does not match')
            end;
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
                self.number_of_clusters = max(self.number_of_clusters, max(self.clusters_labels{subject_id})) ;
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
            for current_subject_index = 1 : length(subject_list)
                current_subject = subject_list(current_subject_index) ;
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
        
        function new_atlas_cluster_label = pre_process_atlas_clusters_labels(self, old_atlas_cluster_label)
            new_number_of_clusters = 0;
            new_atlas_cluster_label = zeros(size(old_atlas_cluster_label)) ;
            vertices_already_relabel = false(size(old_atlas_cluster_label)) ;
            new_atlas_cluster_label(old_atlas_cluster_label == 0) = 0 ;
            vertices_already_relabel(old_atlas_cluster_label == 0) = true ;      
            for current_vertex = 1 : self.number_of_vertices
                if ~vertices_already_relabel(current_vertex)
                    current_vertex_label = old_atlas_cluster_label(current_vertex);
                    new_number_of_clusters = new_number_of_clusters + 1 ;
                    all_searching_vertices = current_vertex ;
                    new_searching_vertices = current_vertex ;
                    while ~ isempty(new_searching_vertices)
                        potential_searching_vertices = [] ;
                        for current_searching_vertex_index = 1 : length(new_searching_vertices)
                            current_searching_vertex = new_searching_vertices(current_searching_vertex_index) ;
                            potential_searching_vertices = union(potential_searching_vertices, self.surface_connectivity_list{current_searching_vertex});
                        end
                        potential_searching_vertices = setdiff(potential_searching_vertices, all_searching_vertices) ;
                        potential_searching_vertices(old_atlas_cluster_label(potential_searching_vertices) ~= current_vertex_label) = [] ;
                        potential_searching_vertices(vertices_already_relabel(potential_searching_vertices)) = [];
                        new_searching_vertices = potential_searching_vertices ;
                        all_searching_vertices = union(all_searching_vertices, new_searching_vertices) ;
                    end
                    if length(all_searching_vertices) >= 50
                        new_atlas_cluster_label(all_searching_vertices) = new_number_of_clusters ;
                        vertices_already_relabel(all_searching_vertices) = true ;
                    else
                        new_number_of_clusters = new_number_of_clusters - 1;
                        new_atlas_cluster_label(all_searching_vertices) = 0 ;
                        vertices_already_relabel(all_searching_vertices) = true ;
                    end
                end
            end
        end
        
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
           for current_subject_index = 1 : length(subject_list)
               current_subject = subject_list(current_subject_index) ;
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
           for current_subject_index = 1 : length(subject_list)
               current_subject = subject_list(current_subject_index) ;
               tmp_average_time_series = zeros(self.number_of_clusters, self.number_of_time_samples) ;
               for current_cluster = 1 : self.number_of_clusters
                   current_cluster_time_series = self.fmri_time_series{current_subject}(self.clusters_vertices_lists{subject_id}{current_cluster}, :) ;
                   tmp_average_time_series(current_cluster, :) = mean(current_cluster_time_series, 1) ;
               end;
               tmp_average_time_series(isnan(tmp_average_time_series)) = 0 ;
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
           for current_subject_index = 1 : length(subject_list)
                current_subject = subject_list(current_subject_index) ;
                tmp_correlation_matrix = cov(self.clusters_average_time_series{current_subject}') ;
                tmp_correlation_matrix(isnan(tmp_correlation_matrix)) = 0 ;
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
            for current_face_index = 1 : length(boundary_faces_local)
                current_face = boundary_faces_local(current_face_index) ;
                boundary_rois_local = union(boundary_rois_local, C(current_face, :)) ;
                boundary_vertices_index(self.cortical_surface.faces(IC == current_face, :)) = true ;
            end
            boundary_vertices_index(self.clusters_labels{subject_id} == 0) = false ;
            boundary_rois_local(boundary_rois_local == 0) = [] ;
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
                a_matrix = diag((diag(l_matrix) + sqrt(diag(l_matrix).^2 + 4 * rho))./ 2./ rho) ;
                tmp_all_inverse_correlation_matrices(:, :, current_subject) = q_matrix * a_matrix * q_matrix' ;
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
            tmp_all_inverse_correlation_matrices(isnan(tmp_all_inverse_correlation_matrices)) = 0 ;
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
                neighbor_rois = union(self.clusters_labels{subject_id}(current_vertex), self.clusters_labels{subject_id}(self.surface_connectivity_list{current_vertex})) ;
                neighbor_rois = intersect(neighbor_rois, self.boundary_rois{subject_id}) ;
                neighbor_rois_index = zeros(size(neighbor_rois)) ;
                for current_roi_index = 1 : length(neighbor_rois)
                    neighbor_rois_index(current_roi_index) = find(self.boundary_rois{subject_id} == neighbor_rois(current_roi_index));
                end;
                neighbor_rois_index(neighbor_rois_index == 0) = [] ;
                for current_roi_index_index = 1 : length(neighbor_rois_index)
                    current_roi_index = neighbor_rois_index(current_roi_index_index) ;
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
                for neighbor_vertex = self.surface_connectivity_list{current_vertex}
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
            rois_clusters = rois_results * (1 : size(rois_results, 2))' ;
            vertices_clusters = vertices_results * (1 : size(vertices_results, 2))' ;
            for current_roi = 1 : number_of_boundary_rois
               candidate_vertices = self.boundary_vertices{subject_id}(vertices_clusters == rois_clusters(current_roi)) ;
               self.clusters_labels{subject_id}(candidate_vertices) = self.boundary_rois{subject_id}(current_roi) ;
            end
        end
        
        function estimate_parcellation(self, subject_id, lambda)
            similarity_matrix = self.compute_similarity_matrix(subject_id, lambda) ;
            graph_cut_results = ncutW(similarity_matrix, length(self.boundary_rois{subject_id})) ;
            self.relabel_boundary_vertices_from_graph_cut_results(subject_id, graph_cut_results) ;
        end;
        
        function subject_instance = create_subject_instance(self, subject_id)
            subject_instance = struct ;
            subject_instance.surface_connectivity_list = self.surface_connectivity_list ;
            subject_instance.boundary_rois = self.boundary_rois{subject_id} ;
            subject_instance.boundary_vertices = self.boundary_vertices{subject_id} ;
            subject_instance.clusters_average_time_series = self.clusters_average_time_series{subject_id} ;
            subject_instance.fmri_time_series = self.fmri_time_series{subject_id} ;
            subject_instance.clusters_inverse_correlation_matrices = self.clusters_inverse_correlation_matrices{subject_id} ;
            subject_instance.clusters_vertices_lists = self.clusters_vertices_lists{subject_id} ;
            subject_instance.clusters_labels = self.clusters_labels{subject_id} ;
            subject_instance.number_of_clusters = self.number_of_clusters ;
        end;
        
        function figure_handle = visualize_both_hemi_cortical_surface(self, subject_id, field_name)
            if ~ exist('field_name', 'var')
                field_name = 'clusters_labels' ;
            end
            if strcmp(field_name, 'clusters_labels')
                face_color_method = 'flat';
            else
                face_color_method = 'interp' ;
            end
            figure_handle = figure ;
            patch('faces', self.cortical_surface.faces, 'vertices', self.cortical_surface.vertices, 'facevertexcdata', self.(field_name){subject_id}, 'facecolor', face_color_method, 'edgecolor', 'none');
            axis equal; axis off; axis vis3d; view(90, 0); camlight; lighting phong; material dull; colormap(jet);
        end
        
        function average_homogeneity = validate_average_homogeneity(self, subject_id)
            assemble_homogeneity = 0;
            number_of_empty_clusters = 0 ;
            for current_cluster = 1 : self.number_of_clusters
                tmp_fmri_time_series = self.reference_data{subject_id}(self.clusters_labels{subject_id} == current_cluster, :) ;
                if ~ isempty(tmp_fmri_time_series)
                    tmp_correlation_matrix = corr(tmp_fmri_time_series') ;
                    tmp_correlation_matrix(isnan(tmp_correlation_matrix)) = 0 ;
                    tmp_correlation_matrix = tmp_correlation_matrix - diag(diag(tmp_correlation_matrix)) ;
                    if size(tmp_correlation_matrix, 1) > 1
                        assemble_homogeneity = assemble_homogeneity + sum(sum(tmp_correlation_matrix)) / ((size(tmp_correlation_matrix, 1) - 1) * (size(tmp_correlation_matrix, 1) - 1)) ;
                    end
                else
                    number_of_empty_clusters = number_of_empty_clusters + 1 ;
                end
            end
            average_homogeneity = assemble_homogeneity / (self.number_of_clusters - number_of_empty_clusters);
        end
        
        function average_variance = validate_average_variance(self, subject_id, task_pair_id)
            max_cluster_mean = 0;
            max_cluster_id = 0;
            for current_cluster = 1 : self.number_of_clusters
                tmp_z_score = self.task_z_scores{subject_id}{task_pair_id}(self.clusters_labels{subject_id} == current_cluster, :);
                if ~ isempty(tmp_z_score)
                   current_cluster_mean = mean(tmp_z_score) ; 
                end
                if current_cluster_mean > max_cluster_mean
                   max_cluster_mean = current_cluster_mean ;
                   max_cluster_id = current_cluster ;
                end
            end
            tmp_z_score = self.task_z_scores{subject_id}{task_pair_id}(self.clusters_labels{subject_id} == max_cluster_id, :) ;
            if ~ isempty(tmp_z_score)
                tmp_z_score_variance = std(tmp_z_score) ;
            end
            average_variance = tmp_z_score_variance ;
        end
        
        function group_connectivity_variance = compute_group_connectivity_variance(self)
            group_connectivity_matrix = zeros(self.number_of_clusters, self.number_of_clusters, self.number_of_subjects) ;
            for current_subject = 1 : self.number_of_subjects
                self.compute_clusters_correlation_matrices(current_subject) ;
                group_connectivity_matrix(:, :, current_subject) = self.clusters_correlation_matrices{current_subject} ;
            end
            assemble_connectivity_variance = sum(sum(var(group_connectivity_matrix, 1, 3))) ;
            group_connectivity_variance = assemble_connectivity_variance / self.number_of_clusters / self.number_of_clusters ;
        end
                    
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
                lambda2 = 0.25 ;
            end;
            %%
            for current_iteration = 1 : max_iterations
                disp(['Now starting the ' num2str(current_iteration) 'th iteration']) ;
                self.history{current_iteration} = self.clusters_labels ;
                self.boundary_vertices = cell(self.number_of_subjects, 1);
                self.boundary_faces = cell(self.number_of_subjects, 1);
                self.boundary_rois = cell(self.number_of_subjects, 1);
                for current_subject = 1 : self.number_of_subjects
                    self.compute_clusters_vertices_lists(current_subject) ;
                    self.compute_clusters_average_time_series(current_subject) ; 
                    self.compute_clusters_correlation_matrices(current_subject) ;
                end
                self.estimate_clusters_inverse_correlation_matrices_group(lambda2, 1) ;
                tmp_subjects_instances = cell(self.number_of_subjects, 1) ;
                tmp_all_clusters_labels = cell(self.number_of_subjects, 1) ;
                for current_subject = 1 :self.number_of_subjects
                    self.find_boundary_vertices(current_subject);
                    tmp_subjects_instances{current_subject} = self.create_subject_instance(current_subject) ;
                end
                disp('Finish estimating group inverse correlation matrices') ;
                parfor current_subject = 1 : self.number_of_subjects
                    display(['Start refining ' num2str(current_subject) 'th subject.']);
                    tmp_all_clusters_labels{current_subject} = ooc_estimate_parcellation(tmp_subjects_instances{current_subject}, lambda1) ;
                end;
                self.clusters_labels = tmp_all_clusters_labels ;
            end;           
        end;
    end;
end
