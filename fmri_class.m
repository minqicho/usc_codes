classdef fmri_class < handle
    
    properties
        cortical_surface ;
        fmri_time_series ;
        clusters_labels ;
        number_of_vertices ;
        number_of_faces ;
        number_of_time_samples ;
        number_of_clusters ;
    end
    
    methods
        
        function self = fmri_class()
            self.number_of_vertices = [] ;
        end
        
        function data_validity = check_data_validity(self, number_of_vertices_local)
            data_validity = false ;
            if isempty(self.number_of_vertices)
                data_validity = true ;
                self.number_of_vertices = number_of_vertices_local ;
            elseif self.number_of_vertices == number_of_vertices_local
                data_validity = ture ;
            end
        end
                
        function set_fmri_time_series(self, fmri_time_series_local)
            [number_of_vertices_local, number_of_time_samples_local] = size(fmri_time_series_local) ;
            if self.check_data_validity(number_of_vertices_local)
                self.fmri_time_series = fmri_time_series_local ;
                self.number_of_time_samples = number_of_time_samples_local ;
            elseif self.check_data_validity(number_of_time_samples_local)
                self.fmri_time_series = fmri_time_series_local' ;
                self.number_of_time_samples = number_of_vertices_local ;
            else
                error('Number of vertices does not match')
            end
        end
        
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
            end
        end
        
        function set_clusters_labels(self, clusters_labels_local, number_of_clusters_local)
            [number_of_vertices_local, second_dimension] = size(clusters_labels_local) ;
            if (self.check_data_validity(number_of_vertices_local)) && (second_dimension == 1)
                self.clusters_labels = clusters_labels_local ;
            elseif (self.check_data_validity(second_dimension)) && (number_of_vertices_local == 1)
                self.clusters_labels = clusters_labels_local ;
            else
                error('Number of vertices does not match')
            end
            if exist('number_of_clusters_local', 'var')
                self.number_of_clusters = number_of_clusters_local ;
            else
                self.number_of_clusters = max(self.clusters_labels) ;
            end
        end 
    end
end
