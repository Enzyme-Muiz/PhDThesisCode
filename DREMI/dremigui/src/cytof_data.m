%How do I extend this to time-series data?
%Need help extending this to time series. I need to make a struct of data
%qualities of time series data is that there arent equally as many cells
%in each time point so the data unless its chopped can't be in one array
%how do I initialize the marker channels?
   
classdef cytof_data
   
    properties    
        
        %basic data 
        data 
        data_backup
        
       
        wells
        well_key
        
        
        header
        %doublets or debris
        %event information 
        debris
       
        event_nums
        invalid 

        beads
        event_gate
        name_channel_map
        channel_name_map

        
        DREMI_adjacency
        DREMI_sig_adjacency
        DREMI_sig_edges
        DREMI_channels
        activity_matrix
        DREMI_noise_threshold
        DREVI_colormode
        
        
        
        ord_clusters
        time_pred_clusters
        
        barcodes
       
        
        
        
        
        
        
        
        bead_nlogl
        bead_post
        gmfit
        bead_gmfit
        dna_gmfit
        %should I remove the doublets? yes probably. 
        
        


        num_time_points
        cell_length_col 
        
        %stims
        %drugs
        %doses
        time_points
        
        %channel information     
        markers
        channels
        proteins
        surface_channels
        internal_channels
        barcoding_channels
        marker_channels
        dna_channels
        bead_channels
        eventnum_channel
        time_col
        
        %computed stats per channel
        means
        corrs
        vars
        skews
        kurtosis
        peaks
        maxes 
        mins
        


        channel_densities
        channel_density_xvals
       
        

        
        %noise and measurement effects 
        channel_decays
        mean_decay
        total_decay
        bc_gmm_means 
        bc_gmm_sigmas
        debris_gmm_means
        debris_gmm_sigmas
        
        
        interpolated_bead_X;
        interpolated_bead_Y;
        
        
        %stuff that used to be in bnorm
        bead_data
        acquisition_times
        acquisition_rate
        acquisition_rate_index
       
        channel_averages
        channel_intercepts
        %if I use the channel averages then I have to compute them. 
        %I should do it on the fly
       
        data_recon
        fluctuation_vector
        dna_low
        
        start_times
        end_times
       
        initial_y
        initial_x
        
        shift_x
        shift_y
        
       
        smooth_decays
        
        
        smooth_corrs
        recon_corrs
        
        smooth_index
        smooth_data
        smooth_diffs  
        smooth_data_recon 
        smooth_fluctuation_vector
        
        
        
        %tsne
        tsne_mapped_data
        sampled_events
        tsne_gated_indices
        

        
    end
    
    methods
        %input output 
        function obj = cytof_data( filename, name2, well_key)
           
            
           obj.DREVI_colormode = 1;
           obj.DREMI_noise_threshold = 0.85;
            if(strcmp(filename,'empty'))
                
                sprintf('no file specified: %s', filename)
                %have to do everything else manually 
            return;
            end
            [obj.data, obj.header, ~] = fca_readfcs(filename); %using arg1
            use_name = 0;
            if (nargin > 1)
                
                if(strcmp(name2,'name'))
                    use_name = 1;
                end
                
            end
                
           
            if(use_name ==0)
                obj = obj.compute_name_channel_maps();    
            else
               obj = obj.compute_name_channel_maps('name'); 
            end
            
            
            [~,c] = size(obj.data);
            obj.marker_channels = linspace(1,c,c);
            
           
            obj.time_col = obj.find_channel_prefix('time'); 
            if(length(obj.time_col)>0)
               
                obj.marker_channels(find(obj.marker_channels == obj.time_col)) = [];  
            end
           
                
            obj.cell_length_col = obj.find_channel_prefix('cell_length'); 
            if(length(obj.cell_length_col)>0)
               
                obj.marker_channels(find(obj.marker_channels == obj.cell_length_col)) = [];
            end
                
           
            obj.dna_channels = obj.find_channel_prefix('dna');
            [~, num_dna]=size(obj.dna_channels);
            for j=1:num_dna
                    
                obj.marker_channels(find(obj.marker_channels == obj.dna_channels(j))) = [];         
                    
            end
            
            obj.barcoding_channels = obj.find_channel_prefix('bc');
            [~, num_barcodes]=size(obj.barcoding_channels)
           
            for j=1:num_barcodes
                    
              obj.marker_channels(find(obj.marker_channels == obj.barcoding_channels(j))) = [];    
                         
            end
            
            obj.bead_channels = obj.find_channel_prefix('bead');
            [~, num_bead]=size(obj.bead_channels);
            for j=1:num_bead
                    
                obj.marker_channels(find(obj.marker_channels == obj.bead_channels(j))) = [];         
                    
            end
          
            trash_channel = obj.find_channel_prefix('barcode');
            [~, num_trash]=size(trash_channel);
            for j=1:num_trash
                    
                obj.marker_channels(find(obj.marker_channels == trash_channel(j))) = [];         
                    
            end
          
             %using arg4
            
            %specifyed barcoding channels (in ORDER)
            if (nargin > 2)
                
                obj.well_key = well_key;
                
            end
            
%             %event number column is already computed and specified 
%             if (nargin > 2)
%                 
%                obj.event_no_col = event_no_col;
%                obj.marker_channels(find(obj.marker_channels == obj.event_no_col)) = [];    
%                
%                
%             end    
%          
            obj.time_points = ones(size(obj.data,1),1);
        end
        
        function channel_indices =  find_channel_prefix(obj, prefix_string)
            
           num_channels = size(obj.data,2); 
           channel_indices = [];
           for i=1:num_channels
               channel_name = obj.channel_name_map(i);
               index = strfind(channel_name,prefix_string);
               if(isempty(index{1}))
                   continue;
               end
               
               if(index{1}(1)~=1)
                   continue;
               end
              channel_indices = [channel_indices i] ;
               
           end
            
        end
        
        function write_data(obj, filename, varargin)
            
            %have to figure out what marker_names and channel names are 
            optargin = size(varargin, 2);
            [num_events, num_markers] = size(obj.data);
            marker_names = cell(num_markers,1);
            channel_names = cell(num_markers,1);
            
            for i=1:num_markers
            
                marker_names(i) = obj.channel_name_map(i);
                channel_names(i) = obj.channel_name_map(i);
                
            end

                
                data = obj.data;
                      fca_writefcs(filename, data, marker_names, channel_names);
    
            
            
        end
        
        
        function obj = compute_name_channel_maps(obj, varargin)
            
            s = sprintf('in name channel map! \n')
        
            use_name = 0;
            optargin = length(varargin);
            for i=1:optargin
               
                if(strcmp('name',varargin{i}))
                    use_name = 1
                end
            end
            
            header_size = length(obj.header.par);
            obj.name_channel_map = containers.Map(); 
            obj.channel_name_map = cell(header_size,1);
            for i=1:header_size 
        
                  hname = obj.header.par(i).name2;
                  if(use_name)
                      hname = obj.header.par(i).name;
                  end
                  
                  if(length(strfind(hname, '('))==0)
                    
                
                      %stuff like the DNA channels do not have the isotope
                      %names
                  
                      if(length(hname)==0)
                           hname = sprintf('channel%d',i);
                            %if its still null for some reason just use the
                            %channel number 
                      end
                     hname = lower(hname);
                     obj.name_channel_map(hname) = i;
                     obj.channel_name_map{i} = hname;
                  
                  else
                       n=strfind(hname, '(');
                
                       if(n==1)
                           %using isotope name 
                        
                            n_end = strfind(hname, ')');
                            real_name = hname(n+1:n_end-1);
                       
                       else
                       
                            real_name = hname(1:n-1);
                       
                       end
                   
            
                        if(length(real_name)==0)
                             real_name = sprintf('channel%d',i);
                             %if its still null for some reason just use the
                             %channel number 
                        end
                        real_name = lower(real_name);
                        obj.name_channel_map(real_name) = i;
                        obj.channel_name_map{i} = real_name;
                 end
              
           end  
     
            
        end
     
        function obj = compute_name_channel_map_from_array(obj, varargin)
            
           
            
            header_size = length(obj.header);
            obj.name_channel_map = containers.Map(); 
            obj.channel_name_map = cell(header_size,1);
            for i=1:header_size 
        
                  hname = obj.header{i};
                  
                  
                  if(length(strfind(hname, '('))==0)
                    
                
                      %stuff like the DNA channels do not have the isotope
                      %names
                  
                      if(length(hname)==0)
                           hname = sprintf('channel%d',i);
                            %if its still null for some reason just use the
                            %channel number 
                      end
                     hname = lower(hname);
                     obj.name_channel_map(hname) = i;
                     obj.channel_name_map{i} = hname;
                  
                  else
                       n=strfind(hname, '(');
                
                       if(n==1)
                           %using isotope name 
                        
                            n_end = strfind(hname, ')');
                            real_name = hname(n+1:n_end-1);
                       
                       else
                       
                            real_name = hname(1:n-1);
                       
                       end
                   
            
                        if(length(real_name)==0)
                             real_name = sprintf('channel%d',i);
                             %if its still null for some reason just use the
                             %channel number 
                        end
                        real_name = lower(real_name);
                        obj.name_channel_map(real_name) = i;
                        obj.channel_name_map{i} = real_name;
                 end
              
           end  
     
            
         end
               
               
       function [name] = get_name_from_channel(obj, channel)
            
           name = obj.channel_name_map(channel); 
            
        end
        
        function [channel] = get_channel_from_name(obj, name)
            
            name = lower(name);
            channel = obj.name_channel_map(name);
        end
        
        function obj = add_data_matrix_header(obj, data, header, wells)
            
           obj.data = data;
           obj.header = header;
           %obj = obj.compute_name_channel_maps();
           obj = obj.compute_name_channel_map_from_array();
            
                [~, c] = size(obj.data);
                for i=1:c
                
                    if i==obj.cell_length_col 
                        continue;
                    end
               
                    if i==obj.time_col 
                        continue;
                    end
                
                    obj.marker_channels = [obj.marker_channels i];
                end
            
            [~, num_dna]=size(obj.dna_channels);
                for j=1:num_dna
                    
                    obj.marker_channels(find(obj.marker_channels == obj.dna_channels(j))) = [];         
                    
                end
            
            [~, num_barcodes]=size(obj.barcoding_channels);
           
                for j=1:num_barcodes
                    
                    obj.marker_channels(find(obj.marker_channels == obj.barcoding_channels(j))) = [];    
                         
                end
                
              
             %obj.marker_channels(find(obj.marker_channels == obj.event_no_col)) = [];    
             obj.wells = wells; 
             %debarcoding info so that I do not have to repeat. 
            
        end
        
        
        function obj = add_dna_channels(obj,dna_channels)
       
            %specifying the DNA channels 
                obj.dna_channels = dna_channels;
                [~, num_dna]=size(dna_channels);
                for j=1:num_dna
                    
                    obj.marker_channels(find(obj.marker_channels == dna_channels(j))) = [];         
                    
                end
                
        end
             
        
        function obj = add_bead_channels(obj,bead_channels)
            
            obj.bead_channels = bead_channels;
            [~, num_beads]=size(bead_channels);
           
%                 for j=1:num_beads
%                     
%                     obj.marker_channels(find(obj.marker_channels == bead_channels(j))) = [];    
%                          
%                 end
            
        end
        
        %specify eventnum channel
        
        function obj = add_eventnum_channel(obj, eventnum_channel)
            
            obj.eventnum_channel = eventnum_channel;
            obj.marker_channels(find(obj.marker_channels == eventnum_channel)) = [];      
                   
        end
        
        function obj = append_cdata_object(obj, cdataobj)
            
            new_data = cdataobj.data;
            obj.data = [obj.data; new_data];
            
            
        end
        
        function obj = add_data(obj, filename, varargin)
            
             %assuming same header for this data otherwise use a higher
             %level construct 
             %ayyo I guess the panels are not the same for these 
             %oh well i guess I'll have to make sure I know what the panels
             %are at some point
             filename
             
             [new_data, ~ , ~] = fca_readfcs(filename); 
                
           
             
             optargin = size(varargin, 2); 
             
             if(optargin == 1)
                 
                 if(strcmp(varargin{1}, 'fix_times'))
                  
                     [old_data_size, ~] = size(obj.data)
                     final_acquisition_time = obj.data(old_data_size,1)
                
                    [new_data_size, ~] = size(new_data)
                    
                    for i=1:new_data_size
                     
                        new_data(i,1) = new_data(i,1)+final_acquisition_time;
                    
                    end
                 
                    obj.data = [obj.data; new_data];
                 end
                 
                 if(strcmp(varargin{1}, 'add timepoints'))
                    
                     
                     
                     obj.data = [obj.data; new_data];
                     
                     time_point_end = length(obj.time_points);
                     new_time_point = obj.time_points(time_point_end)+1;
                     time_points = ones(size(new_data,1),1);
                     time_points = time_points .* new_time_point;
                     obj.time_points = [obj.time_points; time_points];
                    
                     
                 end
             else
                 
                   obj.data = [obj.data; new_data];
                   
             end
             
             
             
        end

        
        

        
         function dna_val = debris_cluster_value(obj,dna_index,eventnum)
          
            dna_channel = obj.dna_channels(dna_index);
           
            prob0 = normpdf(obj.data(eventnum,dna_channel), obj.debris_gmm_means(dna_index,1), obj.debris_gmm_sigmas(dna_index,1));
            prob1 = normpdf(obj.data(eventnum,dna_channel), obj.debris_gmm_means(dna_index,2), obj.debris_gmm_sigmas(dna_index,2));
            
            if(prob0 > prob1) 
                dna_val = 0;
            elseif (prob1 > prob0) 
                dna_val = 1;
            else
                dna_val = -1;
            end
            
        end
        

       function obj = compute_debris_gmm(obj)
        clear X;
        
        [~,num_DNA] = size(obj.dna_channels);
        obj.debris_gmm_means = zeros(num_DNA,2); 
            
            for i=1:num_DNA
               
                DNA_col_data = obj.data(:,obj.dna_channels(i));
                [~, M, V, ~] = EM_GM_fast(DNA_col_data,2); 
                
                if (M(1) < M(2))
                    obj.debris_gmm_means(i,1) = M(1);
                    obj.debris_gmm_means(i,2) = M(2);
                    obj.debris_gmm_sigmas(i,1) = V(1,1,1);
                    obj.debris_gmm_sigmas(i,2) = V(1,1,2);
                else
                    obj.debris_gmm_means(i,1) = M(2);
                    obj.debris_gmm_means(i,2) = M(1);
                    obj.debris_gmm_sigmas(i,1) = V(1,1,2);
                    obj.debris_gmm_sigmas(i,2) = V(1,1,1);
                    
                end    
                
                
            end
        
            
       end
       
       
       
        function obj = compute_stats(obj)
            
           obj = obj.compute_means();
           obj = obj.compute_corrcoefs();
           obj = obj.compute_vars();
           obj = obj.compute_distro_stats();
           obj = obj.compute_maxes_mins();
        end
        
        
        

        

        
        
        function obj = compute_means(obj)
            
               
                [~,num_markers] = size(obj.marker_channels);
                 obj.means = zeros(1,num_markers);
                for i=1:num_markers
                    c=obj.marker_channels(i);
                    obj.means(i) = median(obj.data(:, c));
                 end
        end
        
        function obj = compute_maxes_mins(obj)
            
               
                [~,num_markers] = size(obj.marker_channels);
                 obj.mins = zeros(1,num_markers);
                 obj.maxes = zeros(1,num_markers);
                for i=1:num_markers
                    c=obj.marker_channels(i);
                    obj.mins(i) = min(obj.data(:, c));
                    obj.maxes(i) = max(obj.data(:,c));
                    
                 end
        end
        
        function obj = compute_distro_stats(obj)
            
               
                [~,num_markers] = size(obj.marker_channels);
                 
                obj.skews = zeros(1, num_markers);
                obj.peaks = zeros(1, num_markers);
                obj.kurtosis = zeros(1, num_markers);
                for i=1:num_markers
                    c=obj.marker_channels(i);
                    obj.skews(i) = skewness(obj.data(:, c))
                    obj.kurtosis(i) = kurtosis(obj.data(:,c));
                    [f, xi] = ksdensity(obj.data(:,c));
                    obj.peaks(i)= length(findpeaks(f));
                 end
        end
        
        
        
        
        function channel_mean = get_channel_mean(obj, channel_name)
        
                channel = obj.name_channel_map(channel_name); 
            
                channel_mean = mean(obj.data(:,channel));
            
        end
        
         function channel_var = get_channel_var(obj, channel_name)
        
                channel = obj.name_channel_map(channel_name); 
            
                channel_mean = var(obj.data(:,channel));
            
        end
        
        


        function obj = compute_corrcoefs(obj)
           
            [~,num_markers] = size(obj.marker_channels);
            for i=1:num_markers
                for j = 1: num_markers
                   
                    x = obj.marker_channels(i);
                    y = obj.marker_channels(j);
                    R = corrcoef(obj.data(:, x), obj.data(:, y));
                    obj.corrs(i,j) = R(1,2);
                    
                end
            end
        end
        
        function obj = compute_vars(obj)
            
            [~,num_markers] = size(obj.marker_channels);
            for i=1:num_markers
                c=obj.marker_channels(i);
                x = obj.data(:,c);
                obj.vars(i) = var(x);
                
            end
            
        end

        
        function obj = transform_data(obj,c)
            
            obj.data_backup = obj.data;
            
            am_being_transformed = sprintf('I am being transformed! \n')
           

            
            [~,num_markers] = size(obj.marker_channels);
            for i=1:num_markers
                d=obj.marker_channels(i);
                
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
                
            end
             [~,bc_channels] = size(obj.barcoding_channels);
            for i=1:bc_channels
                d=obj.barcoding_channels(i);
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
              
            end
            
             [~,dn_channels] = size(obj.dna_channels);
            for i=1:dn_channels
                d=obj.dna_channels(i);
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
               
            end
            
        end
        
        function obj = add_computed_channel(obj, channel_name, channel_data)
            
           num_channels = size(obj.data,2);
           obj.data = [obj.data channel_data];
           obj.marker_channels = [obj.marker_channels num_channels+1];
           channel_num = num_channels+1;
           obj.name_channel_map(channel_name) = channel_num;
           obj.channel_name_map{channel_num} = channel_name;
            
        end
            
        function obj = untransform_data(obj)
            
            obj.data = obj.data_backup;
            obj.data_backup = [];
            
            
        end
        
        function obj = remove_low_dna(obj)
        
               obj = obj.compute_dna_gmm_clusters(2, obj.dna_channels);
               obj.debris = obj.dna_low;
               obj = obj.remove_debris();
            
        end
        
        
        
        function obj = clean_beads(obj)
            
           %basically take the beads you get from identify_beads and discard anything with anything with other channels active
          num_events = length(obj.beads);
          for j = 1:num_events 
           
            if(obj.beads(j) == 0)
                continue;
            end
            is_bead = 1;
           
            for i=1:length(obj.marker_channels)
               
                marker_channel = obj.marker_channels(i); 
                
                is_bead_channel = find(obj.bead_channels==marker_channel);
                if(length(is_bead_channel)>0)
                    continue;
                end
              
                marker_channel_value = obj.data(j,marker_channel);
                if(marker_channel_value>0)
                    is_bead = 0;
                    break;
                end
               
               
            end
            
            %do it for dna channels too 
            for i=1:length(obj.dna_channels)
                
                dna_channel = obj.dna_channels(i);
                dna_channel_value = obj.data(j,dna_channel);
                if(dna_channel_value>0)
                    is_bead = 0;
                    break;
                end
                
                
            end
           
            obj.beads(j) = is_bead;
          end
          
          
          
        end
           
            
       
        
        function obj = identify_beads(obj)
           
            
           found_bead_cluster = 0;
           
           
           filter_channels = [obj.bead_channels obj.dna_channels];
           %filter_channels = union(filter_channels, obj.marker_channels)
           
           %dna channels and bead channels, the dna channels should be very
           %low or negative in fact, adding other channels also  
           
           bead_cluster_index = 0;
           
           
            %doublets should be eliminated if I take only cells with 
            %bead-bead 
           sprintf(' in bead identify ')
           
           obj = obj.compute_dna_gmm_clusters(2, obj.dna_channels);
           obj.invalid = ~obj.dna_low;
           
           
           %obj = obj.compute_dna_gmm_clusters(2, 4);
           %obj.invalid = obj.invalid | obj.dna_low
           sprintf('done doing the dna clusters \n')
           
        
           %figure out an approximate percentage of beads automatically 
           
           potential_bead_indices = find(obj.dna_low == 1);
           num_beads = 0;
           %obj.beads = obj.dna_low;
           
           size(potential_bead_indices)
           
           for i=1:length(potential_bead_indices)
              
               index = potential_bead_indices(i);
               is_bead = 1;
               for j=1:length(obj.bead_channels)
                  
                   if(obj.data(index,obj.bead_channels(j))<0)
                       
                      is_bead = 0;
                      break;
                   end
               end
               
               if(is_bead) 
                   
                   num_beads = num_beads + 1;
                   %obj.beads(index) = 1;
               
               end
           end
           [num_events, ~] = size(obj.data); 
           
           approx_percentage = num_beads/ num_events
           approx_num_beads = num_beads
           
           
           for i=2:10 
               
               %precluster the DNA channels and only do it on those. 
               
               
               obj = obj.compute_bead_gmm_clusters(i, filter_channels);
               
               sprintf('done with bead gmm clusters for k= %d \n', i)
               
               if(length(obj.bead_gmfit.mu)==0)
                   %I think that this means it did not converge
                   continue;
               end
               
               %these were sorted this time 
               %so now I have to find the bead channels and the crappy
               %nonbead channels
               bead_indices =zeros(1, length(obj.bead_channels));
              
               for j=1:length(obj.bead_channels)
                   
                   bead_indices(j) = find(filter_channels == obj.bead_channels(j))
                    
               end
               
               nonbead_channels = setdiff(filter_channels, obj.bead_channels)
               nonbead_indices=zeros(1, length(nonbead_channels));
               
               for j=1:length(nonbead_channels)
                   
                   nonbead_indices(j) = find(filter_channels == nonbead_channels(j))
               end
               
               mean_beadsum = sum(obj.bead_gmfit.mu(:,bead_indices),2);
               mean_nonbeadsum = sum(obj.bead_gmfit.mu(:,nonbead_indices),2);
               
               [~,max_bead_index] = max(mean_beadsum);
              
               
               [~,min_nonbead_index] = min(mean_nonbeadsum);
              
               
               if(max_bead_index ~= min_nonbead_index)
                   
                  %sprintf('WARNING max bead and min non beads are not on the same cluster. \n')
                   
               end
               
               bead_cluster_index = max_bead_index;
               %this one is better because the other stuff should mostly be
               %debris
               
                
               cluster_sizes = obj.compute_bc_cluster_sizes(obj.beads)
               [num_events, ~] = size(obj.data);
               percentage = cluster_sizes(bead_cluster_index)/num_events;
               
               if(percentage < approx_percentage)
                   found_bead_cluster = 1;
                   s = sprintf('found a bead cluster! %d components %f percent of events \n', i, percentage)
                   break;
               end
               
           end
            
           if(found_bead_cluster == 0)
               
               s = sprintf('Did not find a bead cluster ! \n')
           else
               
               s = sprintf('Found a bead cluster ! \n')
               num_events = length(obj.beads);     
               
                for i=1:num_events
              
                   if(obj.beads(i) == bead_cluster_index)
                       %just add a criteria about how strong you want that
                       %cluster potential to be 
                       if(obj.bead_post(i,bead_cluster_index)>.7)
                           if(obj.dna_low(i) == 1)
                               obj.beads(i) = 1;
                               sprintf('found a bead');
                           else
                               obj.beads(i) = 0;
                           end
                       else
                           obj.beads(i) = 0;
                       end
                   else
                       obj.beads(i) = 0;
                   end
                end
               

%                 bead_data = obj.data(bead_indices,obj.bead_channels);
%                 all_data = obj.data(:, obj.bead_channels);
%                 within_cluster_distances = mahal(bead_data, bead_data);
%                 distance_threshold = 2*std(within_cluster_distances);
%                 distances_to_beads = mahal(all_data, bead_data);
%                 
%                
%                 
%                 
%                 new_bead_indices = find(distances_to_beads<distance_threshold);
%                 
%                 for m=1:length(new_bead_indices)
%                   
%                     index = new_bead_indices(m);
%                     
%                     if(obj.beads(index)==1)
%                         continue;
%                     end
%                     
%                     if(obj.dna_low(index)==1)
%                         obj.beads(index)=1;
%                     end
%                 end
             
                
           end
           
           obj.clean_beads();
           obj.debris = obj.dna_low;
            
        end
        
        
       
        
        function bead_indices = get_bead_indices(obj)
           
            bead_indices = find(obj.beads==1);
            
        end
        
        function obj = correct_channel_by_marker(obj, channel, windowsize)
            
            
           smooth_data = zeros(1,num_events-window_size);
           smooth_index = zeros(1,num_events-window_size);
            
           for i=1:num_events-window_size
            
               smooth_data(i) = median(event_data(i:i+window_size,channel));
               smooth_index(i) = median(event_data(i:i+window_size,1));
                
           end
            
           
        end 
        
        function obj = plot_bead_channels_vs_dna(obj)
            
                dna_channel = obj.dna_channels(1);
                
           for i=1:length(obj.bead_channels)
              
               subplot(1,5,i)
               bead_channel = obj.bead_channels(i);
               data = [obj.data(:,bead_channel) obj.data(:,dna_channel)];
           
               [bandwidth,density,X,Y]=kde2d(data);
           

                contour(X,Y,density,30), hold on, 
                plot(data(:,1),data(:,2),'b.','MarkerSize',5)
            
                bead_indices = find(obj.beads==1);
            
                data = [obj.data(bead_indices,bead_channel) obj.data(bead_indices,dna_channel)];
           
                [bandwidth,density,X,Y]=kde2d(data);
           

                contour(X,Y,density,30), hold on,
                plot(data(:,1),data(:,2),'r.','MarkerSize',5)
               
               
               
               
           end
            
            
        end
        
        function obj = plot_bead_channel_density(obj, channel1, channel2)
            
          
            
            data = [obj.data(:,channel1) obj.data(:,channel2)];
           
           [bandwidth,density,X,Y]=kde2d(data);
           

            contour(X,Y,density,30), hold on, 
            plot(data(:,1),data(:,2),'b.','MarkerSize',5)
            
            bead_indices = find(obj.beads==1);
            
            data = [obj.data(bead_indices,channel1) obj.data(bead_indices,channel2)];
           
           [bandwidth,density,X,Y]=kde2d(data);
           

            contour(X,Y,density,30), hold on,
            plot(data(:,1),data(:,2),'r.','MarkerSize',5)
            
            
        end
        
        function obj = manual_bead_identify(obj, filter_channel1, filter_channel2)
           
           fhandle = figure;
           
           obj.plot_2d_channel_density(filter_channel1, filter_channel2);
           
           rect = getrect(fhandle)
           left = rect(1);
           bottom = rect(2);
           right = rect(3);
           top = rect(4);
           
           %format for the rectangle is left bottom right top
           
           channel1 = obj.name_channel_map(filter_channel1);
           channel2 = obj.name_channel_map(filter_channel2);
           
           [num_events, ~] = size(obj.data);
           obj.beads = zeros(num_events,1);
           
           for i=1:num_events 
              
              channel1_data = obj.data(i,channel1);
              channel2_data = obj.data(i, channel2);
              if ((left<channel1_data) && (channel1_data<left+right))
                  
                  if ((bottom<channel2_data) && (channel1_data<bottom+top))
                    obj.beads(i) = 1;
                  end
              end
               
           end
           
           num_beads = size(find(obj.beads == 1))
           bead_percentage = num_beads/num_events 
            
        end
        
        function obj = compute_dna_gmm_clusters(obj, num_clusters, filter_channels)
            
            
            clustering_data = obj.data(:,filter_channels);
            options = statset('MaxIter',100);
            
            gmfit = gmdistribution.fit(clustering_data,num_clusters,'options',options);
            
            
            obj.dna_gmfit = gmfit;
           
            [dna_low,nlogl,post] = gmfit.cluster(clustering_data);
          
            
            %wait i don't know which cluster is which 
            %otherwise I may have to switch 0's and ones 
            
            dna_cluster = 0;
            min_sum = 100000;
            
            for i=1:length(gmfit.mu)
                
                if(sum(gmfit.mu(i))<min_sum)
                    
                    min_sum = sum(gmfit.mu(i));
                    dna_cluster = i;
                end
               
            end
            
            dna_cluster
            
            for i=1:length(dna_low)
                
               if(dna_low(i) == dna_cluster)
                   
                  dna_low(i) = 1; 
               else
                   
                  dna_low(i) = 0; 
               end
                
            end
            
            
%             for i=1:length(dna_low)
%                 
%                  if(dna_low(i) == 0)
%                      continue;
%                      
%                  end
%                
%                 if(post(i) < .97)
%                     
%                    
%                     dna_low(i) = 0;
%                     
%                end
%                 
%             end
         
            
            obj.dna_low = dna_low;
            
        end
        
        
          function obj=compute_bead_gmm_clusters(obj, num_clusters, filter_channels)
            
            
            %n-dimensional gaussian mixture model
            % Inputs:
            %   X(n,d) - input data, n=number of observations, d=dimension of variable
            %   k - maximum number of Gaussian components allowed
            %   ltol - percentage of the log likelihood difference between 2 iterations ([] for none)
            %   maxiter - maximum number of iteration allowed ([] for none)
            %   pflag - 1 for plotting GM for 1D or 2D cases only, 0 otherwise ([] for none)
            %   Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)
            obj.beads = [];
            
            bead_data=obj.data(:,filter_channels);
            %gchannels = [3 4];
            %bead_data=obj.data(:,gchannels);
            invalids = find(obj.invalid==1);
            valids = find(obj.invalid==0);
            bead_data(invalids,:)=[];
            
            
            %learning without the obviously invalid ones 
            
            %d = number of barcoding channels here. 
            
            % Ouputs:
            %   W(1,k) - estimated weights of GM
            %   M(d,k) - estimated mean vectors of GM
            %   V(d,d,k) - estimated covariance matrices of GM
            %   L - log likelihood of estimates
            
            options = statset('MaxIter',400);
            
            gmfit = gmdistribution.fit(bead_data,num_clusters,'options',options);
           
            
            obj.bead_gmfit = gmfit;
           
            %bead_data = obj.data(:,filter_channels);
            
            [beads,obj.bead_nlogl,bead_post] = gmfit.cluster(bead_data);
            
            num_valids = length(valids);
            [num_total_events,~] = size(obj.data);
            obj.beads = zeros(num_total_events,1);
            [~,cols] = size(bead_post);
            obj.bead_post = zeros(num_total_events, cols); 
            
            for i=1:num_valids
                obj.beads(valids(i)) = beads(i);
                obj.bead_post(valids(i),:) = bead_post(i,:);
            end
            obj.beads(invalids) = 0; %reset the invalids to 0 
            num_invalids = length(invalids);
            for i=1:num_invalids
                obj.bead_post(i,:) = bead_post(1,:);
            end
            %filling it with some random crap
            
          end
          
          function obj = remove_beads(obj)
             
              bead_indices = find(obj.beads==1);
              
              bead_data = obj.data(bead_indices, obj.bead_channels);
              
              obj.data(bead_indices,:)=[];
              obj.beads(bead_indices)=[];
%               obj = obj.remove_debris();
%               %currently removing all debris including beads
%              
%              
%               distance_from_beads = mahal(obj.data(:,obj.bead_channels), bead_data);
%               
%               distance_from_itself = mahal(bead_data, bead_data);
%               
%               threshold = max(distance_from_itself);
%               
%               bead_doublet_indices = find(distance_from_beads<threshold);
%               
%               obj.data(bead_doublet_indices,:)=[];
              
          end



          

          
          

        function write_data_matrix(obj, filename, channels)
            
            
            M = obj.data(:,channels);
            size_vector = size(M);
            dlmwrite(filename, size_vector, 'delimiter', '\t');
            dlmwrite(filename, M,'delimiter','\t', '-append');
            
            
            
        end
        
        function write_data_matrix_csv(obj, filename)
            
            %Do I want to write out all the channels, I guess so right 
           [num_elements, num_channels] = size(obj.data);
           
            fid = fopen(filename,'w');
            for i=1:num_channels-1
           
               
                fprintf(fid,'%s,',obj.channel_name_map{i});
       
               
           
            end 
            
            fprintf(fid,'%s \n',obj.channel_name_map{num_channels});
            fclose(fid);
            
           for i=1:num_channels
           
            
            dlmwrite(filename, obj.data,'delimiter',',', '-append');
           
           end 
            
            
        end

        
        

        function obj = compute_decay(obj)
            
            [~,n] = size(obj.data);
            for i=1:n
                channel = i;
                p = polyfit(obj.data(:, obj.eventnum_channel), obj.data(:,channel),1);
                obj.channel_decays(i) = p(1);
            end
            
            obj.mean_decay = mean(obj.channel_decays);
            obj.total_decay = sum(obj.channel_decays);
        end
        
        function P = compute_channel_decay(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            
            P = polyfit(obj.data(:, obj.eventnum_channel), obj.data(:,channel),1);
            
          
        end
        
        
        function [smooth_index, smooth_data] =  plot_channel_smooth(obj, channel_name, window_size, varargin)
            
             channel = obj.name_channel_map(channel_name);
            
            [num_events, n] = size(obj.data);
            smooth_data = zeros(1,num_events-window_size);
            for i=1:num_events-window_size
                smooth_data(i) = median(obj.data(i:i+window_size,channel));
            end
            smooth_index=linspace(1,num_events-window_size,num_events-window_size);
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(smooth_index,smooth_data,varargin{1});
            else
                plot(smooth_index,smooth_data);
            end
            
            
        end
        
        function [smooth_index, smooth_data] =  plot_channel_vs_time_smooth(obj, channel_name, window_size, varargin)
            
            
            channel = obj.name_channel_map(channel_name);
            [num_events, n] = size(obj.data);
            smooth_data = zeros(1,num_events-window_size);
            smooth_index = zeros(1,num_events-window_size);
            for i=1:num_events-window_size
                smooth_data(i) = median(obj.data(i:i+window_size,channel));
                smooth_index(i) = median(obj.data(i:i+window_size,1));
            end
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(smooth_index,smooth_data,varargin{1});
            else
                plot(smooth_index,smooth_data);
            end
            
            
        end
           

        
        function obj = remove_nonbeads(obj)
            
            non_bead_events = find(obj.beads==0);
            bead_events = find(obj.beads==1);
            
            sprintf('num bead events %d num nonbead events %d ', length(bead_events), length(non_bead_events))
            
            obj.data(non_bead_events,:) = [];
            
        end
        
        function obj = remove_nongated(obj)
            
            non_gated_events = find(obj.event_gate==0);
            gated_events = find(obj.event_gate==1);
            
            sprintf('num gated events %d num nongated events %d ', length(gated_events), length(non_gated_events))
            
            obj.data(non_gated_events,:) = [];
            
        end
        
        function obj = remove_debris(obj)
            
            
            %the low dna clusters that are not beads
            debris_indices = find(obj.debris==1);
            obj.data(debris_indices,:) = [];
        end
        

        function plot_channel_means(obj, channel_name, window_size)
            
            channel = obj.name_channel_map(channel_name);
            [num_events, ~] = size(obj.data);
            data_means=zeros(num_events-window_size,1);
            for i=1:num_events-window_size
                
                data_means(i) = mean(obj.data(i:i+window_size,channel));
            end
            
            plot(data_means);
            
        end
        
        function obj = compute_debris(obj)
        
            %For now lets try removing events with -DNA and -length
            [num_events, ~] = size(obj.data); 
            
            obj.debris = zeros(1, num_events); 
            for i=1:num_events
            
                debris_cluster1 = obj.debris_cluster_value(1,i);
                debris_cluster2 = obj.debris_cluster_value(2,i);
                
                if (debris_cluster1 == 0) || (debris_cluster2 == 0)
                    obj.debris(i) = 1;
                else
                    obj.debris(i) = 0;
                end
            end
            
                %If it seemed like it was generated from the first rather
                %than the second cluster then discarded as the debris 
            
        end
        
        function plot_2d_channel_density(obj, channel1_name, channel2_name, varargin)
            
           channel1 = obj.name_channel_map(channel1_name);
           channel2 = obj.name_channel_map(channel2_name); 
           
           data = [obj.data(:,channel1) obj.data(:,channel2)];
           limits = [];
           for i=1:length(varargin)
               
                
               if(strcmp(varargin{i}, 'limits'))
                  
                   limits = varargin{i+1};
                   
               end
           end
           
           
           %[bandwidth,density,X,Y]=kde2d(data,256, [min(obj.data(:,channel1)) min(obj.data(:,channel2))], [max(obj.data(:,channel1)) max(obj.data(:,channel2))]);
           maxx = max(obj.data(:,channel1));
          
           if(length(limits)>0)
                
               [bandwidth,density,X,Y]=kde2d(data,256, [0 0], [limits(3) limits(4)]);
           else
               [bandwidth,density,X,Y]=kde2d(data,256, [0 0], [max(obj.data(:,channel1)) max(obj.data(:,channel2))]);
           end
           
           %http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
           %Dani recommends this method of density estimation 
           %http://code.google.com/p/danapeerlab/source/browse/trunk/freecell/src/axes.py
           %python for doing this.
           %Dani says that the contour3 stuff shoudl work. 
           
           %leaving out optional arguments
           %[bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
            
           slices_size = maxx/8;
           
         
           
            %surf(X,Y,density) 
            %view([0,70])
            %colormap hot 
            %hold on 
            %alpha(.8)
            %set(gca, 'color', 'blue');
            %plot(data(:,1),data(:,2),'w.','MarkerSize',5)
           
            optargin = length(varargin);
            if(optargin == 0)
            
                plot(data(:,1),data(:,2),'b.','MarkerSize',5) 
                hold on, 
                contour(X,Y,density,30), 
            end
            %
            %plot(data(:,1),data(:,2),'b.','MarkerSize',5)
            
            for i=1:length(varargin)
               
               
               if(strcmp(varargin{i}, 'imagesc'))
                   
%                   density_filtered = density>.15;
%                   density = density.*density_filtered;
%                   j = jet;
%                   j(1,:) = [ 1 1 1 ];
%                   colormap(j);
                  colormap(jet);
                  imagesc(X(1,:),Y(:,1),density); 
                  set(gca,'YDir','normal');
                   set(gca,'XTick',[]);
                   set(gca,'YTick',[]);
                  %plot_as_vertical_lines([1 2 3 4 5 6 7].*slices_size, 'w');
               end  
               if(strcmp(varargin{i}, 'contour'))
                   
                   plot(data(:,1),data(:,2),'b.','MarkerSize',5) 
                   hold on, 
                   contour(X,Y,density,30), 
                   set(gca,'XLim',[0 max(data(:,1))]);
                   set(gca,'YLim',[0 max(data(:,2))]);
               end
                
            end
           
        end
        
        
        
        
        function obj =  tsne_map_data(obj,  varargin)
        
            optargin = length(varargin);
            
            num_events = size(obj.data,1);
           
            
              channels = [];
            if(optargin == 0)
                channels = obj.marker_channels;
            else
               
              
                channel_names = varargin{1};
                for i=1:length(channel_names)
                   
                    channels = [channels obj.name_channel_map(channel_names{i})];
                    
                end
%                channels = varargin{1};
                
            end
            channels
            obj.tsne_mapped_data = fast_tsne(obj.data(:,channels));
                
           
           
        end
        
        function draw_tsne(obj,varargin)
            
            
               if(length(varargin)>0)
                
                %coloring_channel = varargin{1};
                %channel = obj.name_channel_map(coloring_channel);
                color_channel_data=varargin{1}; 
                
                %color_channel_data = obj.data(:,channel);
                scatter(obj.tsne_mapped_data(:,1),obj.tsne_mapped_data(:,2),14 ,color_channel_data,'fill');
               else
                   
                scatter(obj.tsne_mapped_data(:,1),obj.tsne_mapped_data(:,2),14);
            
               
               end
               
               
               
               
%             
%         
%             [~, density, x, y] = kde2d([obj.tsne_mapped_data(:,1) obj.tsne_mapped_data(:,2)], 256);
% 
%             contour(x, y, density, 5);
%             
%             hold
%             
%             if(length(varargin)>0)
%                 
% %                coloring_channel = varargin{1};
% %                channel = obj.name_channel_map(coloring_channel);
%                
% %                color_channel_data = obj.data(:,channel);
%                 color_channel_data = varargin{1};
%                
%                rx = linspace( 0, 1 ); % bins used in color coding
%                colors = jet( length( rx ) - 1 ); % coloring scheme
%                color_channel_data = color_channel_data/max(color_channel_data);
% 
% 
% 
%             	for i = 1:length(rx)-1
%                     rx_indices = find( rx(i) <= color_channel_data & color_channel_data <= rx( i + 1 ) );
%               		plot( obj.tsne_mapped_data( rx_indices, 1 ), obj.tsne_mapped_data( rx_indices, 2 ), '.', 'Color', colors( i, : ) );
%                     
%                 end
% 	
%             end
               
        end
      
        
        
        
        function compute_and_graph_render_edges_dremi(obj, edges, nodes)
           
            
           adjMatrix = zeros(length(nodes), length(nodes));
           edges_indexed = zeros(length(edges),2);
           for i = 1:length(edges)
            
               node_index1 = find(strcmp(edges{i,1},nodes));
               node_index2 = find(strcmp(edges{i,2},nodes));
               adjMatrix(node_index1, node_index2) = 1;
               edges_indexed(i,1) = node_index1;
               edges_indexed(i,2) = node_index2;
               
                
           end 
           gObj = biograph(adjMatrix,nodes); 
          
           set(gObj,'ShowArrows','on');
           set(gObj, 'ShowWeights','on');
           set(gObj,'EdgeFontSize',12);
           
           for i=1:length(edges) 
            
                
                dremi_values(i) = obj.compute_dremi(edges{i,1},edges{i,2},.80);
            
           end
           
           range = .65;
           load 'continuous_BuPu9.mat';
           colormap(continuous_BuPu9); 
           
          
           
           for i = 1:length(dremi_values)
               
               [color_value] = get_color_value(dremi_values(i), .65, continuous_BuPu9);
               set(gObj.edges(i),'LineColor',color_value);
               set(gObj.edges(i),'LineWidth',2.0);
               set(gObj.edges(i),'Weight',dremi_values(i));
           end
           
           all_nodes = 1:length(nodes); 
           set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
           set(gObj.nodes(all_nodes),'Color',[1 1 1]);
           set(gObj.nodes(all_nodes),'Shape','ellipse');
           set(gObj.nodes(all_nodes),'LineWidth',1.1);
           set(gObj.nodes(all_nodes),'fontSize',10);
           
           view(gObj);
           
        end
        
        
        
        function [channel_names] = get_marker_channel_names(obj)
            
            channel_names = cell(1,length(obj.marker_channels));
            
            for i=1:length(channel_names)
                
               channel_names{i} = obj.channel_name_map{obj.marker_channels(i)}; 
            end
            
        end
        
        function obj = tsne_gate(obj, num_rects)

            fhandle = figure;
            mapped_data = obj.tsne_mapped_data;
           
            [~, density, x, y] = kde2d([mapped_data(:,1) mapped_data(:,2)], 256);

            contour(x, y, density, 12);
    
           
            gated_indices = [];
            num_events = size(mapped_data,1);
    
    
            for j=1:num_rects      
        
             rect = getrect(fhandle)
             left = rect(1);
            bottom = rect(2);
            width = rect(3);
            height = rect(4);
    
                for i=1:num_events
    
                    if((mapped_data(i,1)>left)&&(mapped_data(i,1)<left+width))
                        if((mapped_data(i,2)>bottom)&&(mapped_data(i,2)<bottom+height))
                
                        gated_indices = union(gated_indices ,[i]);
                
                        end
                    end
        
                end
        
            end
            obj.tsne_gated_indices = gated_indices
           
            percentage = length(obj.tsne_gated_indices)/length(obj.sampled_events)
            
        end
         
        function visualize_markers(obj,varargin)
            
            marker_channels = obj.marker_channels;
            
            optargin = length(varargin);
            
            if(optargin>0)
               
                draw_channel_names = varargin{1};
                
                for i=1:length(draw_channel_names)
                    
                   new_channel = obj.name_channel_map(draw_channel_names(i))
                   marker_channels = [marker_channels new_channel];
                    
                end
                
            end
            
            ncols = 5;
            nrows = ceil(length(marker_channels)/ncols);

            for i=1:length(marker_channels)
        
                channel_name = obj.channel_name_map{marker_channels(i)};
                subplot(nrows,ncols,i);   
                obj.plot_channel_density(channel_name,'b');
        
                xlabel(channel_name, 'fontsize', 11);
       
      
        
            end
             
        end
        
        
    function [signaling_direction_vector]  = visualize_cluster_markers(obj, cluster_indices, varargin)
            
          num_events = size(obj.data, 1);
           
          non_cluster_indices = linspace(1, num_events, num_events);
          non_cluster_indices(cluster_indices) = [];
        
          signaling_direction_vector = [];

          optargin = length(varargin);
          channels = [];
          if(optargin == 0)
             
              channels = obj.marker_channels;
          else
              
              channel_names = varargin{1};
             for i=1:length(channel_names)
                
                 channel_no = obj.name_channel_map(channel_names{i});
                 channels = [channels channel_no];
                 
             end
              
              
          end
          
          ncols = 5;
          nrows = ceil(length(channels)/5);
          
          for i=1:length(channels)
        
            subplot(nrows,ncols,i);   
            channel = channels(i);
      
            channel_name = obj.channel_name_map{channel}; 
            marker_data = obj.data(:,channel);
    
            non_cluster_data = marker_data(non_cluster_indices);
            cluster_data = marker_data(cluster_indices);
       
            %signaling_direction_vector{1,i} = channel_name;
            population_median = median(obj.data(:,channel));
            %signaling_direction_vector{2,i} = median(cluster_data)-median(population_median)/median(population_median);
            signaling_direction_vector(i) = median(cluster_data)-median(population_median)/median(population_median);
       
            [f, xi] = ksdensity(non_cluster_data);
            plot(xi,f,'LineWidth',1.3);
            hold 
            f2 = ksdensity(cluster_data,xi);
            plot(xi,f2, 'r','LineWidth',1.3);
            xlabel(channel_name, 'fontsize', 11);
       
      
        end
    
            
            
        end
        
        function plot_2d_channel_scatter(obj, channel1_name, channel2_name,varargin)
            
           channel1 = obj.name_channel_map(channel1_name);
           channel2 = obj.name_channel_map(channel2_name); 
           channel1_data = obj.data(:,channel1);
           channel2_data = obj.data(:,channel2);
           low_x_indices = find(channel1_data<0);
           channel1_data(low_x_indices)=[];
           channel2_data(low_x_indices)=[];
           
           low_x_indices = find(channel2_data<0);
           channel1_data(low_x_indices)=[];
           channel2_data(low_x_indices)=[];
           
           optargin = length(varargin);
           if(optargin>0)
                 plot(channel1_data,channel2_data,varargin{1},'*');
           else
           
               plot(channel1_data, channel2_data,'*');
           end
           minx = min(obj.data(:,channel1));
           miny = min(obj.data(:,channel2));
           maxx = max(obj.data(:,channel1));
           maxy = max(obj.data(:,channel2));
           xlim([0 maxx]);
           ylim([0 maxy]);
           
           set(gca,'XTick',[]);
           set(gca,'YTick',[]);
           %box on 
           
        end
        
        
        function pairwise_visual_with_density(obj, channel1_name, channel2_name, varargin)
            
            subplot(2,1,1);
            [points_x, points_y, normalized_density] = obj.pairwise_visualize(channel1_name, channel2_name, varargin);
            
            channel1 = obj.name_channel_map(channel1_name);
            
            f = ksdensity(obj.data(:,channel1), points_x);
            
            maxf = max(f);
            f = f/maxf;
            
            subplot(2,1,2);
            
            plot(points_x, f);
            
            
        end



      
        
        


     function [obj] = prune_transitive_edges(obj)
            %this is a really limited case of dismissing edges
        
        
       
       
        for i=1:length(obj.DREMI_channels)
            for j=1:length(obj.DREMI_channels)
            
             
              if(obj.DREMI_sig_adjacency(i,j)==0)
                  continue;
              end
              
              %found an edge
              
              for k=1:length(obj.DREMI_channels)
              
                  if(obj.DREMI_sig_adjacency(j,k)==0)
                      continue;
                  end
                      
                  if(obj.DREMI_sig_adjacency(i,k)==0)
                      continue;
                  end
                  
                  %found a potential case of dismissable transitivity
                  
                  if(obj.DREMI_sig_adjacency(i,k)<obj.DREMI_sig_adjacency(j,k))
                           
                     obj.DREMI_sig_adjacency(i,k) = 0;
                    
                  end
           
              
              
              end
            end
            
            
        end
        
        
        
        
        end
        
        
        
        function draw_denovo_DREMI_graph(obj)
            
          
           binary_adjacency = (obj.DREMI_sig_adjacency>0);
           colsum = sum(binary_adjacency,1);
           rowsum = sum(binary_adjacency,2);
           index = 0;
           node_map = zeros(1,length(colsum));
           
           nodes_to_keep=[];
           nodes_to_delete=[];
           for i=1:length(colsum)
               if((colsum(i)>0)|(rowsum(i)>0))
                   index = index+1;
                   node_map(i)=index;
                   nodes_to_keep = [nodes_to_keep i];
               else
                   nodes_to_delete = [nodes_to_delete i];
               end
           end
           
           reduced_adjacency = binary_adjacency(nodes_to_keep, nodes_to_keep);
           num_nodes = length(reduced_adjacency);
           reduced_channels = obj.DREMI_channels(nodes_to_keep);
           reduced_sig_adjacency = obj.DREMI_sig_adjacency(nodes_to_keep, nodes_to_keep);
           
           gObj = biograph(reduced_adjacency,reduced_channels); 
           
           range = .7;
           load 'continuous_BuPu9.mat';
           colormap(continuous_BuPu9); 
           all_nodes = 1:num_nodes; 
           set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
           set(gObj.nodes(all_nodes),'Color',[1 1 1]);
           set(gObj.nodes(all_nodes),'Shape','ellipse');
           set(gObj.nodes(all_nodes),'LineWidth',1.1);
           set(gObj.nodes(all_nodes),'fontSize',10);
           set(gObj,'ShowArrows','on');
           set(gObj, 'ShowWeights','on');
           set(gObj,'EdgeFontSize',12);
%            all_nodes = 1:length(obj.DREMI_channels); 
%            set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
%            set(gObj.nodes(all_nodes),'Color',[1 1 1]);
%            set(gObj.nodes(all_nodes),'Shape','ellipse');
%            set(gObj.nodes(all_nodes),'LineWidth',1.1);
%            set(gObj.nodes(all_nodes),'fontSize',10);
           
           
           for i = 1:size(reduced_adjacency,1)
               for j = 1:size(reduced_adjacency,2)
                   
                  
                   if(reduced_adjacency(i,j)==0)
                       continue;
                   end
                  
                   
                   edgeID = getedgesbynodeid(gObj,get(gObj.Nodes([i j]),'ID'));
                  
                 [color_value] = get_color_value(reduced_sig_adjacency(i,j), .7, continuous_BuPu9);
                 set(edgeID,'LineColor',color_value);
                 set(edgeID,'LineWidth',2.0);
                 labelstring = sprintf('%.2f',reduced_sig_adjacency(i,j));
                 set(edgeID,'Label',labelstring);
                 set(edgeID,'Weight', reduced_sig_adjacency(i,j));
                 
               end
           end
           
           dolayout(gObj);
           view(gObj);
           
           
        end
       
        
        
        
        
        function [obj, sig_edges, sig_mis, sig_edge_matrix] = write_mi_graph(obj, channel_names, threshold, varargin)
       
               mi_matrix = obj.DREMI_adjacency;
               num_edges = length(find(mi_matrix>threshold)); 
               sig_edges = {num_edges, 2};
               sig_mis = zeros(num_edges,1);
               %file = fopen(filename,'w');
               undirected = 0; 
               choose_direction = 0;
               sig_edge_matrix = zeros(length(channel_names), length(channel_names));
               
               if(~isempty(varargin))
                     if(strcmp(varargin{1},'undirected'))
                            undirected = 1; 
                     end
                     if(strcmp(varargin{1},'choose'))
                            choose_direction = 1; 
                     end
               end
     
              
                current = 1;
                for j=1:length(channel_names)
                    for k=1:length(channel_names)
                        
                        if(undirected ==1)
                            if(j==k)
                                break;
                            end
                        end
                        
                        val = mi_matrix(j,k);
                        if(undirected == 1)
                            val = max(mi_matrix(j,k), mi_matrix(k,j));
                        end    
                        
                        if(val>threshold)
                            
                             if(choose_direction ==1)
                                 if(val<mi_matrix(k,j)) %will only enter one direction 
                                     continue;
                                 end
                             end
                            
                              %fprintf(file,'%s %s %f\n', channel_names{j}, channel_names{k}, val);
                             
                              sig_edges{current,1} = channel_names{j};
                              sig_edges{current,2} = channel_names{k};
                              sig_mis(current) = val;
                              current = current+1;
                              sig_edge_matrix(j,k) = val;
                        end
                        
                    end
                end
                obj.DREMI_sig_adjacency = sig_edge_matrix;
                obj.DREMI_channels = channel_names;
                obj.DREMI_sig_edges = sig_edges;
                %obj.DREMI_sig_edge_mi = sig_mis;
                
        end
        
        function [molecule_edges] = get_edges_for_y_molecule(obj, channel_name)
           
           edge_indices = []; 
           num_edges = length(obj.DREMI_sig_edges);
           for i=1:num_edges
               
                edge = obj.DREMI_sig_edges(i,:);
                if(strcmp(edge{2}, channel_name))
                    edge_indices = [edge_indices i];
                end
               
           end
            
            molecule_edges = obj.DREMI_sig_edges(edge_indices,:);
        end
        
        
        function [auc, eval_points_x, eval_points_y] = compute_edge_auc(obj, channel1_name, channel2_name, varargin)
            
            
            
            eval_points_x = [];
            
            for i=1:length(varargin)
            
                         if(strcmp(varargin{i}, 'eval_points'))
            
                            
                                eval_points_x = varargin{i+1};
                    
            
                         end
                    
                         
            end
            
            
           
            
           [points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
            
            if(length(eval_points_x)==0)
               
                eval_points_x = points_x;
                
            end
              
            eval_points_y = interp1(points_x,points_y,eval_points_x, 'linear', 'extrap');
            
            %area under curve 
            auc = sum(eval_points_y .* abs(eval_points_x(2) - eval_points_x(1)));
            
            %normalize the area under the curve. 
            %auc 
            xrange = max(eval_points_x)-min(eval_points_x);
            yrange = max(eval_points_y)-min(eval_points_y);
            auc = auc/(xrange*yrange);
            
            
        end
            
        function [obj, activity_matrix] = pairwise_auc_compute(obj, channel_names, varargin)
           
            activity_matrix = zeros(length(channel_names), length(channel_names));
            ctrl_specified = 0;
            for i=1:length(varargin)
            
                         if(strcmp(varargin{i}, 'ctrl_data'))
            
                            ctrl_specified = 1; 
                            ctrl_data = varargin{i+1};
                    
            
                         end   
         
           end
            
            
            for i=1:length(channel_names)
                for j=1:length(channel_names)
                    
                    if i==j
                        continue;
                    end
                    
                   
                    if(ctrl_specified==1)
                        
                        %sprintf('computing activity matrices, control specified')
                        
                        [~, points_x1, points_y1] = obj.compute_edge_auc(channel_names{i},channel_names{j});
                        [~, points_x2, points_y2] = ctrl_data.compute_edge_auc(channel_names{i},channel_names{j});
                        
                       if(max(points_x1) < max(points_x2))
                           points_x = points_x1;
                       else
                           points_x = points_x2;
                       end
                        %i dont think union or intersection works. 
                        
                        [auc1, points_x, points_y1] = obj.compute_edge_auc(channel_names{i},channel_names{j}, 'eval_points', points_x);
                        [auc2, points_x, points_y2] = ctrl_data.compute_edge_auc(channel_names{i},channel_names{j}, 'eval_points', points_x);
                        
                         points_y = abs(points_y1-points_y2);
%                         auc = sum(points_y);
                        
                        auc = sum(points_y .* abs(points_x(2) - points_x(1)));
                        xrange = max(points_x)-min(points_x);
                        yrange = max(points_y)-min(points_y);
                        auc = auc/(xrange*yrange);
                        dremi1 = obj.compute_dremi(channel_names{i}, channel_names{j}, 0.8);
                        dremi2 = ctrl_data.compute_dremi(channel_names{i}, channel_names{j}, 0.8);
                        dremi = max(dremi1, dremi2);
                        
                        if(dremi>.15)
                            activity_matrix(i,j) = auc;
                        end
                           
                    else
                        
                        activity_matrix(i,j) = obj.compute_edge_auc(channel_names{i},channel_names{j});
                        
                        
                        
                    end
                    
                end
            end
            
                obj.activity_matrix = activity_matrix;
                sprintf('set object activity matrix to activity matrix')
                
            
            
        end
        
        
        function [obj, mi_matrix ] = pairwise_mi_compute(obj, channel_names, noise_threshold, varargin)
    
           
            ctrl_specified = -1; 
            ctrl_data = [];
            
            obj.DREMI_noise_threshold  = noise_threshold;
            
            mi_matrix = zeros(length(channel_names), length(channel_names));
            
            skiplist = [];
           
            
            for i=1:length(varargin)
            
                         if(strcmp(varargin{i}, 'ctrl_data'))
            
                            ctrl_specified = 1; 
                            ctrl_data = varargin{i+1};
                    
            
                         end
                         if(strcmp(varargin{i}, 'skip_list'))
            
                           
                            skiplist = varargin{i+1};
                    
            
                         end
            
            
         
           end
            
             
            
            for i=1:length(channel_names)
                for j=1:length(channel_names)
                    if(i==j) 
                        continue;
                    end
                    channel1_name = channel_names{i}
                    channel2_name = channel_names{j}
                    
                    skip=0;
                    for k=1:size(skiplist,1)
                   
                        if(strcmp(skiplist(k,1), channel1_name))
                            if(strcmp(skiplist(k,2),channel2_name))
                            
                                skip = 1;
                                break;
                            end
                        end
                    end
                
                    if(skip == 1)
                        continue;
                        %skip this computation
                    end
                
            
                    if(ctrl_specified<0)
                        [mi_matrix(i,j),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold);
                    else
                        [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(obj, channel1_name, channel2_name, 25, 255);
                        [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(ctrl_data, channel1_name, channel2_name, 25, 255);
                        maxy = max(maxy1,maxy2);
                        [mi_matrix(i,j),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold, 'maxy', maxy);
                    end
                    
                    
                    
                    
                
                end
            end
        
            obj.DREMI_adjacency = mi_matrix;
        
            for i=1:length(varargin)
            
                if(strcmp(varargin{i}, 'plot'))
            
          
                    colormap(jet)
                    CLIM = [0 1];
                    imagesc(mi_matrix);
                    set(gca,'ytick',1:length(channel_names));
                    set(gca,'yticklabel',channel_names);
                    xticklabel_rotate([1:length(channel_names)],45,channel_names);
                    colorbar
            
                end
            
            
            
         
            end
        
        
        end
        
        
        
        
        function [obj] = subsample_data(obj, percent_subsampled)
            
            data_length = size(obj.data,1);
            new_sample = randsample(data_length,floor(data_length*percent_subsampled));
            obj.data = obj.data(new_sample,:);
            
            
        end
        
        
        
        function [dremi, pvalue, samples_x, samples_y] = compute_dremi(obj, channel1_name, channel2_name, noise_threshold , varargin)
            
            
             compute_pvalue = 0;
             set_maxy  = 0;
             num_permutations = 0; 
             max_yval = 0; 
             num_slices = 8; 
             non_kde_style = 0;
             min_yval = 0;
             
             for i=1:length(varargin)
                 
                if(strcmp(varargin{i}, 'compute_pvalue'))
                    compute_pvalue = 1;
                    num_permutations = varargin{i+1};
                end
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                   
                    maxy_val = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'num_slices'))
                    num_slices = varargin{i+1};
                end
                if(strcmp(varargin{i}, 'no_kde'))
                   non_kde_style = 1;  
                end
                 
                
             end
             
             use_min_y = max([0, min(obj.data(:, obj.name_channel_map(channel2_name)))]);
             
             if(non_kde_style==0)

                if(set_maxy==0); 
                 [points_x, points_y, point_weights, ~, normalized_density, ~,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', noise_threshold,'no_plot');
                else
                 
                 [points_x, points_y, point_weights, ~, normalized_density, ~,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', noise_threshold,'no_plot', 'MinMaxY', use_min_y, maxy_val);
                end
                
                  total_slice_samples = sum(point_weights) * 1000;
                  samples_x=[];
                  samples_y=[];
             
                  for i = 1:length(points_x)
                 
                      num_repeats = floor(point_weights(i) * 1000);
                      new_samples_x = ones(num_repeats,1).*points_x(i);
                      new_samples_y = ones(num_repeats,1).*points_y(i);
        
                      samples_x = [samples_x; new_samples_x];
                      samples_y = [samples_y; new_samples_y];
                  end
                
             
             
                 data = [samples_x samples_y];
        
%             [points_x, points_y, point_weights, ~, normalized_density, xaxis,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
%             data = [transpose(points_x) transpose(points_y)];
            
            
            %[dremi, entropy_y, cond_entropy_y] = dremi_resampled(data, 8, 8);
               minx = min(data(:,1));
               miny = min(data(:,2));
               maxx = max(data(:,1));
               maxy = max(data(:,2));
               
             else
                 
                [minx, miny, maxx, maxy] = find_data_cutoffs(obj, channel1_name, channel2_name, 50, 255);
                 
                 
             end
            
          
            
            if(set_maxy==1)
                maxy = maxy_val;
            end
            
            if(non_kde_style == 0)
                
                dremi = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices,num_slices);
                
            else
                
                dremi = delta_entropyreweight(obj, channel1_name, channel2_name, num_slices, num_slices, minx, miny, maxx, maxy)
   
            end
            
            pvalue = 0;
            if(num_permutations == 0)
                compute_pvalue = 0; 
            end
            
            if(compute_pvalue ==1)
                s= sprintf('computing pvalue!\n')
                dremi_random = zeros(1,num_permutations);
                total_perms_performed = num_permutations;
                num_greater = 0; 
                for i=1:num_permutations
                    
                    if(non_kde_style==0)
                        
                         dremi_random(i) = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices, num_slices, 'permute_y');
                        
                         
                    else
                        
                        dremi_random(i) = delta_entropyreweight(obj, channel1_name, channel2_name, num_slices, num_slices, minx, miny, maxx, maxy);
                        
                    end
   
                    
                    random_dremi = dremi_random(i)
                    if(dremi_random(i)>dremi)
                        s=sprintf('dremi random greater than dremi')
                        
                        num_greater = num_greater+1;
                        if(num_greater==10)
                            total_perms_performed = i;
                            break;
                        end
                        
                    end
        
                end
                
              
                
                above_dremi = find(dremi_random(1:total_perms_performed)>dremi);
                if(total_perms_performed<num_permutations)
                    pvalue = 10/total_perms_performed;
                else
                    pvalue = (length(above_dremi)+1)/(num_permutations+1);
                    
                end
            end
            
            
   
            
            
        end
    
        
      
        

        
        function [F_fitted, MSE, x, y_fit] = double_sigmoid_fit_edge(obj, channel1_name, channel2_name, varargin)
            
            [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.95, 'no_plot');
%            [x,y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
%             channel1 = obj.name_channel_map(channel1_name);
%                 channel2 = obj.name_channel_map(channel2_name);
%                 x = transpose(obj.data(:,channel1));
%                 y =  transpose(obj.data(:,channel2));
            hold on;
            
            for i=1:length(varargin)
                if(strcmp(varargin{i},'cutoff_data'))
                    
                    x_cutoff = varargin{i+1};
                    index = find(x>x_cutoff,1,'first');
                    x(index:end) = [];
                    y(index:end) = [];
                    
                end
            end
            
            lower_intercept = 0;
            mid_intercept = 1; 
            upper_intercept = 2;
            
            fix_intercepts = 0; 
            for i=1:length(varargin)
                if(strcmp(varargin{i},'fix_intercepts'))
                    
                    fix_intercepts = 1;
                   
                    lower_intercept = varargin{i+1};
                    mid_intercept = varargin{i+2};
                    upper_intercept = varargin{i+3};
                    %slope = varargin{4};
                end
            end
            
            if(fix_intercepts==1)

                   f = @(p,x) (((lower_intercept - mid_intercept) ./ (1 + (x/p(1)).^p(2)))+mid_intercept) + (((mid_intercept - upper_intercept) ./ (1 + (x/p(3)).^p(4)))+upper_intercept);
                   [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1 1 1 1]);
                   
            else
            
              
                f = @(p,x) (((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4))+(((p(4) - p(7)) ./ (1 + (x/p(6)).^p(5)))+p(7));
                %p1 = lower, p4 = mid, p7 = upper p3 = inflection1 p2 =
                %slope 1 p6 = infection p5 = slope 2 
                
               
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[lower_intercept 1 1 mid_intercept 1 1 upper_intercept]);
                
                
               
               
            end
            
            disp(['F = ',num2str(F_fitted)])
            % Plot the data and fit
            %figure(1)
            y_fit = f(F_fitted,x);
            plot(x,y,'*',x,f(F_fitted,x),'r');
            %plot(x,y_fit,'w','LineWidth',2.0);
            
            
        end
        
        
         function [X, Y_vals] = conditional_mean_edge(obj, channel1_name, channel2_name, varargin)
             disp('doing this');
             hold on;
            Limits = [];
            for i=1:length(varargin)
               if(strcmp(varargin{i},'Limits')) 
                    Limits = varargin{i+1};
               end
            end
           
            if(length(Limits)==0)
    
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name);
                 
            else
                
                 [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'Limits', Limits);
                
            end
            
             [X, Y_vals, Y_errs] = plot_dremi_curves(x(:), y(:), 'num_locs', 256, 'normalize', false, 'smooth', 0.38, 'avg_type', 'make_plot', 'gaussian', 'svGolay', false);
             plot(X,Y_vals,'w','LineWidth',3.0);
            %plot(x, smooth(y, 20), 'w', 'LineWidth', 3.0);
            hold off;
         end
        
        function [F_fitted, MSE, x,y_fit, Kd ] = sigmoid_fit_edge(obj, channel1_name, channel2_name, varargin)

            hold on;
            Limits = [];
            for i=1:length(varargin)
               if(strcmp(varargin{i},'Limits')) 
                    Limits = varargin{i+1};
               end
            end
        
            if(length(Limits)==0)
    
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.9, 'no_plot');
                 %[x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'no_plot');
                 channel1 = obj.name_channel_map(channel1_name);
                 channel2 = obj.name_channel_map(channel2_name);
                 x = transpose(obj.data(:,channel1));
                 y =  transpose(obj.data(:,channel2));
            else
                
                 [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts', .9,  'Limits', Limits);
                
            end
           
           
            
            for i=1:length(varargin)
                if(strcmp(varargin{i},'cutoff_data'))
                    
                    x_cutoff = varargin{i+1};
                    index = find(x>x_cutoff,1,'first');
                    x(index:end) = [];
                    y(index:end) = [];
                    
                end
            end
            

            
            %f = @(F,x) y = F(1)./(1+exp(-F(2)*(x-F(3))));
            %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)).^p(5))+p(4);
         
            %p1 = lower asymptote
            %p2 = hill slope
            %p3 = disassociatin constant
            %p4 = upper asymptote
            ysmooth = zeros(1,length(x));
            ysmooth(1) = (y(1)+y(2))/2;
            for i = 2:length(x)-1
                
               ysmooth(i) = (y(i-1)+y(i)+y(i+1))/3;
                
            end
            ysmooth(length(x)) = (y(length(x))+y(length(x)-1))/2;
            
            smooth_diff = diff(ysmooth);
            min_diff = min(smooth_diff);
            std_diff = std(smooth_diff)*2; 
            
            index_first = find(smooth_diff>min_diff+std_diff,1, 'first');
            
            index_last = find(fliplr(smooth_diff)>min_diff+std_diff,1,'first');
            index_last = length(x) - index_last; 
            lower_intercept = mean(y(1:index_first))
            upper_intercept = mean(y(index_last:length(x)))
            x_index_first = x(index_first);
            x_index_last = x(index_last);
            
            R = corrcoef(x,y);
            
            if(R<0)
                temp = upper_intercept; 
                upper_intercept = lower_intercept;
                lower_intercept = temp;
            end
            
 
            fix_intercepts = 0; 
            for i=1:length(varargin)
                if(strcmp(varargin{i},'fix_intercepts'))
                    fix_intercepts = 1
                    lower_intercept = varargin{i+1}
                    upper_intercept = varargin{i+2}
                    %slope = varargin{4};
                end
            end
            
            if(fix_intercepts==1)

                   f = @(p,x) ((lower_intercept - upper_intercept) ./ (1 + (x/p(1)).^p(2)))+upper_intercept;
                   [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1 1]);
                   %f = @(p,x) ((lower_intercept - upper_intercept) ./ (1 + (x/p(1)).^slope))+upper_intercept;
                   %[F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1]);
                   Kd = F_fitted(1)^F_fitted(2);
                   Hc = F_fitted(2);
                   inflection = F_fitted(1);
                   hill_slope = F_fitted(2);
            else
            
                %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4);
                %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x.*p(3)).^p(2)))+p(4);
                f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4);
                %f = @(p,x) ((x.^p(2) ./ (p(3)^p(2) + x.^p(2))).* p(1))+p(4);
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[lower_intercept 1 1 upper_intercept]);
                Kd = F_fitted(3) ^ F_fitted(2);
                Hc = F_fitted(2);
                inflection = F_fitted(3);
                hill_slope = F_fitted(2);
            end
            %slope = F_fitted(2)
            %Ka = F_fitted(3)
            %Kd = F_fitted(3)^F_fitted(2)
            % Display fitted coefficients
            %disp(['F = ',num2str(F_fitted)])
            % Plot the data and fit
            %figure(1)
            y_fit = f(F_fitted,x);
            %plot(x,y,'*',x,y_fit,'r');
            
            plot(x,y_fit,'w','LineWidth',3.0);
            xlabel_string = sprintf('inflection: %.2f hill slope: %.2f ', inflection, hill_slope);
            xlabel(xlabel_string);
            %area(x,y_fit,'FaceColor',[0.65,0.65,0.65])
            %set(gca,'XLim',[Limits(1) Limits(3)]);
            %set(gca,'YLim',[Limits(2) Limits(4)]);
            
            %legend('data','fit')

            hold off;
        end
        

        function [Pfit,Delta, points_x, y_fit] = linear_fit_edge(obj, channel1_name, channel2_name, varargin)
            
            hold on;
            
            Limits = [];
            for i=1:length(varargin)
               if(strcmp(varargin{i},'Limits')) 
                    Limits = varargin{i+1};
               end
            end
            
            if(length(Limits)==0)
    
%                [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.90, 'no_plot');
%                  channel1 = obj.name_channel_map(channel1_name);
%                  channel2 = obj.name_channel_map(channel2_name);
%                  points_x = transpose(obj.data(:,channel1));
%                  points_y =  transpose(obj.data(:,channel2));
                 [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
                
            else
                
                 [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.9, 'Limits', Limits);
                
            end
           
          
            
            %[points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
            
            %need to search for a fit           
           
            [Pfit, S1] = polyfit(points_x,points_y,1);
            [y_fit, Delta] = polyval(Pfit, points_x, S1); 
            %hold on;
       
            %y_fit = f(F_fitted,x);
            %plot(points_x,points_y,'*',points_x,y_fit,'r');
            plot(points_x,y_fit,'w','LineWidth',3.0);
            slope = Pfit(1);
            intercept = Pfit(2);
%             area(points_x,y_fit,'FaceColor',[0.65,0.65,0.65]);
            set(gca,'XLim',[Limits(1) Limits(3)]);
            set(gca,'YLim',[Limits(2) Limits(4)]);
            %[r2 rmse] = rsquare(points_y,y_fit) 
            xlabel_string = sprintf('Slope: %.2f Intercept: %.2f', slope, intercept);
            xlabel(xlabel_string);
            hold off;
            
        end
        
        function [R ] = corrcoef_edge(obj, channel1_name, channel2_name)
            
             %[points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
             channel1 = obj.name_channel_map(channel1_name);
             channel2 = obj.name_channel_map(channel2_name);
             points_x = obj.data(:,channel1);
             points_y = obj.data(:,channel2); 
             
             R = corrcoef(points_x, points_y);
        end
       
        function [split, delta, Y1, Y2, Pfit1, Pfit2] = two_line_fit_edge(obj, channel1_name, channel2_name)
        
            [points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
            
            %need to search for a fit 
            delta = 1000;
            split = 0;
            Pfit1 = [];
            Pfit2 = [];
            Y1 = [];
            Y2 = [];
    
    
    
            num_points = length(points_x);
    
            for i=2:num_points-2
        
           
           
                [P1, S1] = polyfit(points_x(1:i),points_y(1:i),1);
                [P2, S2] = polyfit(points_x(i+1:num_points),points_y(i+1:num_points),1);
            
                [x_eval,Delta1] = polyval(P1, points_x(1:i), S1);
                [y_eval,Delta2] = polyval(P2, points_x(i+1:num_points),S2);
                d1 = mean(Delta1);
                d2 = mean(Delta2);
           
                if(max(d1,d2)<delta)
               
                    delta = max(d1,d2);
                    split = i;
                    Pfit1 = P1;
                    Pfit2 = P2;
                    Y1 = x_eval;
                    Y2 = y_eval;
                end
           
           
            end
    
            hold on;
            %plot(points_x, points_y , '*');
            plot(points_x(1:split),Y1,'w','LineWidth',2.0);
            plot(points_x(split+1:num_points),Y2,'w','LineWidth',2.0);

        end
        
      function [points_x, points_y, point_weights, density, normalized_density, xaxis,yaxis] = pairwise_visualize(obj, channel1_name, channel2_name, varargin)
    
        channel1 = obj.name_channel_map(channel1_name);
        channel2 = obj.name_channel_map(channel2_name);
        X = obj.data(:,channel1);
        Y = obj.data(:,channel2);
        total_cells = size(obj.data,1);  
          
        %done = 0;
            
        num_slices = 256;
        minxval = 0;
        minyval = 0;
        cutoff = 50;
        draw_contour = 0;
        show_density = 0;
        draw_plot = 1;
        avg_pts = 1;
        avg_pts_threshold = .9;
        fix_limits = 0;
        maxyval = max(Y);
        maxxval = max(X);
        fixy = 0;
        visual_threshold = 0; 
        axes_specified = 0; 
        axes_handle = 0;
        
        for i=1:length(varargin)-1
            
            if(strcmp(varargin{i},'axes'))
                
                axes_specified = 1;
                axes = varargin{i+1};
            end
                
            if(strcmp(varargin{i},'Slices'))
                num_slices = varargin{i+1};
            end
            if(strcmp(varargin{i},'MinMaxY'))
               
               minyval = varargin{i+1};
               maxyval = varargin{i+2}; 
               
               fixy = 1;
                
            end
            if(strcmp(varargin{i},'MinMaxX'))
               
               minxval = varargin{i+1}; 
               maxxval = varargin{i+2}; 
               fixx = 1;
                
            end
            if(strcmp(varargin{i},'MinMaxFitX'))
                
               minxval = min(X); 
            end
            if(strcmp(varargin{i},'MinMaxFitY'))
                
               minyval = min(Y); 
            end    
                
            if(strcmp(varargin{i},'Cutoff'))
                
               cutoff = varargin{i+1};
               
            end
            
            if(strcmp(varargin{i},'Minval'))
                
               minxval = varargin{i+1};
               minyval = minxval;
            end
            if(strcmp(varargin{i},'Limits'))
                
                fix_limits = 1;
                limitvector = varargin{i+1};
                minx = limitvector(1);
                miny = limitvector(2);
                maxx = limitvector(3);
                maxy = limitvector(4);
                
                
            end
            if(strcmp(varargin{i},'non_averaged_pts'))
                
               avg_pts = 0;
               avg_pts_threshold = varargin{i+1};
               
            end
            if(strcmp(varargin{i}, 'visual_threshold'))
                visual_threshold = varargin{i+1};
            end
            
            
        end
        
       
        
        for i=1:length(varargin)
            
             if(strcmp(varargin{i},'draw_contour'))
                draw_contour = 1;
             end
             if(strcmp(varargin{i},'show_density'))
                 show_density = 1;
             end
             if(strcmp(varargin{i},'no_plot'))
                 draw_plot = 0;
             end
            
        end
        
     
        if (fix_limits == 0)
                Xmin = minxval;
                Xmax = maxxval;
                Ymin = minyval;
                Ymax = maxyval;
        else
                Xmin = minx;
                Xmax = maxx;
                Ymin = miny;
                Ymax = maxy;
        end
      
        
        if fixy == 1
           Ymin = minyval;
           Ymax = maxyval;
        end
%         if done == 0;
%             R = Xmax - Xmin;
%             
%             dx = R/(num_slices - 1);
%             xmesh1 = Xmin + [0:dx:R];
%             [init_X, ind_hist] = histc(X, xmesh1);
%             z_init_X = (init_X - mean(init_X))/std(init_X);
%             ind_high = find(z_init_X > 3);
%             ind_to_remove = [];
%             if ~isempty(ind_high)
%                 hist_sorted = sort(init_X, 'descend');
%                 for i = 1:length(ind_high)
%                     to_add = ceil(mean(init_X) + mean(init_X)/5); %hist_sorted(length(ind_high) + 1);
%                     ind_orig = find(ind_hist == ind_high(i));
%                     ind_to_remove = [ind_to_remove; randsample(ind_orig, length(ind_orig) - to_add)];
%                 end
%                 X(ind_to_remove) = [];
%                 Y(ind_to_remove) = [];
%             end
%             
%             
%             
%             %Xmin = 0;
%             %Xmax = max(X);
%             %Ymin = 0;
%             %Ymax = max(Y);
%             Ry = Ymax - Ymin;
%             dy = Ry/(num_slices - 1);
%             ymesh1 = Ymin + [0:dy:Ry];
%             [init_Y, ind_hist_y] = histc(Y, ymesh1);
%             z_init_Y = (init_Y - mean(init_Y))/std(init_Y);
%             ind_high_y = find(z_init_Y > 3);
%             ind_to_remove_y = [];
%             if ~isempty(ind_high_y)
%                 hist_sorted_y = sort(init_Y, 'descend');
%                 for i = 1:length(ind_high_y)
%                     %to_add_y = hist_sorted_y(length(ind_high_y) + 1);
%                     to_add_y = ceil(mean(init_Y) + mean(init_Y)/5);
%                     ind_orig_y = find(ind_hist_y == ind_high_y(i));
%                     %ind_to_remove_y = [ind_to_remove_y; randsample(ind_orig_y, ceil(mean(init_X) + 1.5*std(init_X)))];
%                     ind_to_remove_y = [ind_to_remove_y; randsample(ind_orig_y, length(ind_orig_y) - to_add_y)];
%                 end
%                 X(ind_to_remove_y) = [];
%                 Y(ind_to_remove_y) = [];
%             end
%             
%             done = 1;
%       end
     
    
        if(fix_limits == 0)
          
            [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[Xmin Ymin],[Xmax Ymax]);
            
        else
            
            [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minx miny],[maxx maxy]);
            
        end
        
        %rounding down small values 
   
        for i=1:num_slices
            for j=1:num_slices
           
                if(density(i,j)<0.00001) 
                    density(i,j) = 0;
                end
            end
        end
   
        
        %weeding out sparse ends 
        
   
        
          row=0;
   
          total_cells_so_far = 0; 
          while(total_cells_so_far<cutoff)
      
    
            row = row+1;
            total_cells_so_far = length(find(Y>Grid_Y(num_slices-row,1)));
        
        
          end
 
          maxy = Grid_Y(num_slices-row,1);
  
      
         total_cells_so_far = 0;
         start_row = 0;
         while(total_cells_so_far<cutoff)
      
       
       
            start_row = start_row + 1;
            total_cells_so_far = length(find(Y<Grid_Y(start_row,1)));
        
         end
       
 
        miny = Grid_Y(start_row,1);
       
        %row = 0;
        %start_row = 1;
      
       total_cells_so_far = 0;
       col=0;
       while(total_cells_so_far<cutoff)
      
       
       
            col = col+1;
            total_cells_so_far = length(find(X>Grid_X(1,num_slices-col)));
       
       end
   
       maxx = Grid_X(1,num_slices-col);
   
   
      total_cells_so_far = 0;
      start_col=0;
      while(total_cells_so_far<cutoff)
      
       
       
         start_col = start_col+1;
         total_cells_so_far = length(find(X<Grid_X(1,start_col)));
       
      end
      
      minx = Grid_X(1,start_col);
   
  
   
   if(fix_limits == 1)
       
          start_row = 1;
          start_col = 1;
          row = 0;
          col = 0;
   end
      
   if(fixy==1)
       
      start_row = 1; 
      row = 0;
       
   end
   density = density(start_row:num_slices-row,start_col:num_slices-col);
   num_cols = size(density,2);
   num_rows = size(density,1);
   xaxis = Grid_X(1,start_col:num_slices-col);
   yaxis = Grid_Y(start_row:num_slices-row,1);
   
   normalized_density = zeros(num_rows,num_cols);
   prob_normalized_density = zeros(num_rows,num_cols);
   %normalized by column for plotting the data 
   for i=1:num_cols
      
      %normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
      normalized_density(:,i) = density(:,i)/max(density(:,i));
      prob_normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
      
   end
   
   
   
   %now create the side bars 
   
   
   colsum = sum(density,1);
   normalized_colsum = colsum./max(colsum);
 
   rowsum = sum(density,2);
   normalized_rowsum = rowsum./max(rowsum);
   
   
   
   
  %the corner is a fudge 
   
  
  %blueval = min(normalized_colsum);
  blueval = 0;
  corner = ones(11,11).*blueval;
 
  %make the top bar
 
  %yaxis_increment = abs(yaxis(2)-yaxis(1,1));
  yaxis_increment = .01;
  yaxis_top_bar = [];
  top_bar = [];
  zero_vector = zeros(1,length(normalized_colsum));
  for i=1:1
      top_bar = [top_bar; zero_vector]; 
      yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];    
      
  end
  for i=1:10
           top_bar = [top_bar; normalized_colsum];
           yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];    
  end
  
   
  %make the side bar
   %xaxis_increment = abs(xaxis(2)-xaxis(1));
   xaxis_increment = .01;
   xaxis_side_bar = [];
   side_bar = [];
   zero_vector = zeros(length(normalized_rowsum),1);
   
   for i=1:1
      side_bar = [side_bar zero_vector]; 
      xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)];    
      
   end
   
   for i=1:10
       side_bar = [side_bar normalized_rowsum];
       xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)]; 
   end
   
   
   
   %find the trace through the peak regions for the return value
   points_x = [];
   points_y = [];
   point_weights = [];
   if(avg_pts==1)
       for i=1:num_cols
       
        
           
           max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
           % points_y = [points_y mean(Grid_Y(max_indices,i))];
           
           points_x = [points_x xaxis(i)];

           %new_point_y = dot(Grid_Y(start_row+max_indices,start_col+i),normalized_density(max_indices,i));
           points_y = [points_y mean(yaxis(max_indices))];
           %points_y = [points_y new_point_y];
      
       end
       point_weights = ones(1,length(points_y));
   else
       
       for i=1:num_cols
       
           %instead of referring to the grid maybe just take all the points
           %in the high density squares ??
           
          
           
           max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
           % points_y = [points_y mean(Grid_Y(max_indices,i))];
           new_points = ones(1,length(max_indices)).*xaxis(i);
           new_point_weights = transpose(normalized_density(max_indices,i));
           new_point_weights = new_point_weights ./ (sum(new_point_weights));
           points_x = [points_x new_points];

           %points_y(i) = dot(Grid_Y(start_row:255-row,i),orig_normalized_density(:,i));
           y_indices = max_indices;
           new_points_y = transpose(yaxis(y_indices));
           points_y = [points_y new_points_y];
           point_weights = [point_weights new_point_weights];
      
       end
       
   end
   
   
   
   smoothed_normalized_density = zeros(num_rows, num_cols);
   
   for i=1:num_rows
       for j = 2:num_cols-1
            smoothed_normalized_density(i,j) = (normalized_density(i,j-1)+normalized_density(i,j)+normalized_density(i,j+1))/3; 
       end
       
   end
   %imagesc(flipud(Grid_X(:,1)), flipud(transpose(Grid_Y(1,:))), normalized_density);
   
   
   if(visual_threshold>0)
       smoothed_normalized_density = (smoothed_normalized_density>visual_threshold).*smoothed_normalized_density;

   end
   matrix_to_plot = [smoothed_normalized_density side_bar];
   top_bar = [top_bar corner];
   matrix_to_plot = [matrix_to_plot; top_bar];
   
   xaxis_to_plot = [xaxis xaxis_side_bar];
   yaxis_to_plot = [yaxis; yaxis_top_bar];
   
   
   
   
   if(draw_plot)
       obj.DREVI_colormode
       if(obj.DREVI_colormode==1)

            colormap(jet);
        
            if(axes_specified==1)
                 %    imagesc(xaxis, yaxis, smoothed_normalized_density);
                imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot, 'Parent', axes_handle);
                set(axes_handle,'YDir','normal');
               
                %set(axes_handle,'Xtick',1:ceil(max(xaxis)));
                %set(axes_handle,'YTick',1:ceil(max(yaxis)));
                xlabel(channel1_name, 'FontSize', 16, 'FontWeight', 'bold')
                ylabel(channel2_name, 'FontSize', 16, 'FontWeight', 'bold')
            else
            
                imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                set(gca,'YDir','normal');
                
                %set(gca,'XTick',1:ceil(max(xaxis)));
                %set(gca,'YTick',1:ceil(max(yaxis)));
                xlabel(channel1_name, 'FontSize', 16, 'FontWeight', 'bold')
                ylabel(channel2_name, 'FontSize', 16, 'FontWeight', 'bold')
        
            end
            
       end
       
        if(obj.DREVI_colormode==2)
              density_filtered = matrix_to_plot>.6;
              matrix_to_plot = matrix_to_plot.*density_filtered;
             j = jet;
             j(1,:) = [ 1 1 1 ];
            colormap(j);
   
  
      
          if(axes_specified==1)
                %         density_filtered = smoothed_normalized_density>.6;
                %         smoothed_normalized_density = smoothed_normalized_density.*density_filtered;

                %    imagesc(xaxis, yaxis, smoothed_normalized_density);
                %xaxis
                imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot, 'Parent', axes_handle);
                set(axes_handle,'YDir','normal');
                % set(axes_handle,'Xtick',1:ceil(max(xaxis)));
                %set(axes_handle,'YTick',1:ceil(max(yaxis)));
                xlabel(channel1_name, 'FontSize', 12, 'FontWeight', 'bold')
                ylabel(channel2_name, 'FontSize', 12, 'FontWeight', 'bold')
          else
                %xaxis
                imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                set(gca,'YDir','normal');
                % set(gca,'XTick',1:ceil(max(xaxis)));
                %set(gca,'YTick',1:ceil(max(yaxis)));
                xlabel(channel1_name, 'FontSize', 16, 'FontWeight', 'bold')
                ylabel(channel2_name, 'FontSize', 16, 'FontWeight', 'bold')
        
            end
        end  
    %    set(axes_handle,'YDir','normal');
        %set(gca,'XTick',[]);
        %set(gca,'YTick',[]);
        %set(gca, 'XTickLabel','');
        %set(gca, 'YTickLabel','');
       
       
%         set(gca, 'FontSize',16);
%         xlabel(channel1_name);
%         ylabel(channel2_name);
        
     %   hold
   end
   
   if(draw_contour)
       [bandwidth,rdensity,rGrid_X,rGrid_Y]=kde2d([X Y],num_slices+1,[minx miny],[maxx maxy]);
  
       contour(rGrid_X, rGrid_Y, rdensity, 12);
       
   end
   
   if(show_density)
       
       f = ksdensity(X, points_x);
       plot(points_x,f, 'w', 'LineWidth',1.3);
       
   end
   
  
  
   
   for i=1:length(varargin)-1
   
    if(strcmp(varargin{i},'Title'))
   
      
       title(varargin{i+1});
      
   
    end
   end
  
 
      end

      function cluster_data_objects = split_to_clusters(obj, cluster_indices)
        
        num_indices = max(cluster_indices);
        cluster_data_objects = cell(1,num_indices);
        
        for i=1:num_indices
            cluster_data_objects{i} = obj; 
            new_data = obj.data(cluster_indices,:);
            cluster_data_objects{i}.data = new_data;
            
            
        end
          
          
      end
      
      
        function obj = cluster_gate_events(obj, cluster_indices)
            
            new_data = obj.data(cluster_indices,:);
            obj.data = new_data;
        end
        
        function [obj, gated_indices] = threshold_gate_events(obj, channel_name, thresh, greater)
            
            channel = obj.name_channel_map(channel_name);
           
            channel_data = obj.data(:,channel);
            
            if(strcmp(greater,'gt'))
                gated_indices = find(channel_data>thresh);
            end
            if(strcmp(greater,'lt'))
                gated_indices = find(channel_data<thresh);
            end
            
               
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
        end
        
        function [obj, gated_indices] = high_low_gate_events(obj, channel_name, percentage_start, percentage_end)
           
            channel1 = obj.name_channel_map(channel_name);
            
            
            channel_data = obj.data(:,channel1);
            
          
            
            num_percentage_start = floor(percentage_start*length(channel_data))+1;
            num_percentage_end = floor(percentage_end*length(channel_data));
            
                
          
           [~, sorted_indices] = sort(channel_data,'ascend');
                
      
                
            
           gated_indices = sorted_indices(num_percentage_start:num_percentage_end);
           new_data = obj.data(gated_indices,:);
           obj.data = new_data;
            
            
            
            
        end
        
      function [obj, gated_indices] = high_low_gate_events_interval(obj, channel_name, interval_start, interval_end)
           
           channel1 = obj.name_channel_map(channel_name);
            
            
           channel_data = obj.data(:,channel1);
                
           high_indices = find(channel_data>interval_start);
           low_indices = find(channel_data(high_indices)<interval_end);
           gated_indices = high_indices(low_indices);
           
           new_data = obj.data(gated_indices,:);
           obj.data = new_data;
            
            
            
            
      end     
        

        function obj = manual_gate_events_box(obj, channel_name1, channel_name2, left, bottom, width, height)
            
            channel1 = obj.name_channel_map(channel_name1);
            channel2 = obj.name_channel_map(channel_name2);
            channel_data = [obj.data(:,channel1) obj.data(:,channel2)]; 
            candidate_indices = find((channel_data(:,1)>left)&(channel_data(:,1)<left+width));
                
            shorter_data_set = [channel_data(candidate_indices,1) channel_data(candidate_indices,2)];
                
            nested_indices = find((shorter_data_set(:,2)>bottom)&(shorter_data_set(:,2)<bottom+height));
                
            gated_indices = candidate_indices(nested_indices);
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            if(length(obj.time_points)>0)
                new_time_points = obj.time_points(gated_indices,:);
                obj.time_points = new_time_points;
            end
        end
        
        function obj = manual_gate_events(obj, channel_name1, channel_name2, num_rects)
            
            fhandle = figure;
            channel1 = obj.name_channel_map(channel_name1);
            channel2 = obj.name_channel_map(channel_name2);
            channel_data = [obj.data(:,channel1) obj.data(:,channel2)]; 
             [~, density, x, y] = kde2d([channel_data(:,1) channel_data(:,2)], 256);

            plot(channel_data(:,1),channel_data(:,2),'b.','MarkerSize',5)
            hold on,
            contour(x, y, density, 30);
            
            
            gated_indices = [];
            num_events = size(obj.data,1);
        
    
            for j=1:num_rects      
        
                rect = getrect(fhandle)
                left = rect(1)
                bottom = rect(2)
                width = rect(3)
                height = rect(4)
    
                candidate_indices = find((channel_data(:,1)>left)&(channel_data(:,1)<left+width) & (channel_data(:,2)>bottom)&(channel_data(:,2)<bottom+height));
               
                
                if(length(gated_indices)==0)
                    gated_indices = candidate_indices;
                else
                    gated_indices = union(gated_indices,candidate_indices);
                end
                    
        
            end
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
        end
        
        function [F, XI] = plot_channel_density(obj, channel_name, varargin)
        
            XI_given = 0; 
            XI = [];
            for i=1:length(varargin)
                if(strcmp(varargin{i},'XI'))
                    XI_given = 1;
                    XI = varargin{i+1};
                end
            end
            
            channel = obj.name_channel_map(channel_name);
            
            if(XI_given==0)
                [F, XI] = ksdensity(obj.data(:,channel));
            else
                [F] = ksdensity(obj.data(:,channel), XI);
                
            end
            
%             optargin = size(varargin,2);
%             if(optargin==1)
%                 plot(XI,F,varargin{1},'LineWidth',5);
%             else
                 plot(XI,F,'LineWidth',5);
%             end
          
%             data = obj.data(:,channel);
%             sorted_data = sort(data);
%             num_points = length(data);
%             low_index = floor(num_points*.01);
%             high_index = floor(num_points-(num_points*.01));
%             low_value = sorted_data(low_index);
%             high_value = sorted_data(high_index);
%             plot_as_vertical_lines([low_value high_value],'r');
            
%             box on;
%             set(gca,'XTick', []);
%             set(gca, 'YTick', []);
            
        end
        
		function plot_channel_hist(obj, channel_name, num_bins, varargin)
											 
			channel = obj.name_channel_map(channel_name);
			data = obj.data(:,channel);
            optargin = size(varargin,2);
            h = findobj(gca,'Type','patch');
            
            if(optargin==1)
                set(h,'FaceColor',varargin{1},'EdgeColor','w');
                hist(data, num_bins);
            else
                set(h,'FaceColor', 'b','EdgeColor','w');
                hist(data, num_bins);
            end
			
		end
											 
        function [gmfit, cluster_ids, cluster_percents] = compute_channel_gm_mixture(obj, channel_name, K)
            
            channel = obj.name_channel_map(channel_name);
            gmfit = gmdistribution.fit(obj.data(:,channel),K);
            cluster_ids = cluster(gmfit,obj.data(:,channel));
            cluster_percents = zeros(K,1);
            total = size(obj.data,1);
            for i=1:K
               
                cluster_percents(i) = length(find(cluster_ids==i))/total;
            end
           
            %bar(cluster_sizes);
            
        end

        
        
        function [channel_data] = get_channel_data(obj, channel_name)

           channel = obj.name_channel_map(channel_name);
           channel_data = obj.data(:, channel);
           
        end
        





        
        function obj = perform_bead_gate(obj, gate_channels, bead_threshold, nonbead_threshold)
            
            [num_events,~] = size(obj.data)
            bead_threshold
            nonbead_threshold
            num_markers = length(obj.marker_channels);
            num_dna = length(obj.dna_channels);
            obj.beads = ones(1,num_events);
            
 
             %repeat the same thing for dna channels
             
              for i=1:num_events
                 
                 for j=1:num_dna
                    
                     dna_channel = obj.dna_channels(j);
                     if(obj.data(i,dna_channel)>nonbead_threshold) 
                          obj.beads(i) = 0;
                          break;
                     end
                     
                 end
             
              end

           
          %bead channels can overlap with marker channels    
          for i=1:num_events 
              
                 
                 %skip the ones that are already considered as non beads
                 if(obj.beads(i) == 0)
                     continue;
                 end
                 
                 for j=1:length(gate_channels)
                     
                     bead_channel = gate_channels(j);
                     if(obj.data(i,bead_channel)<bead_threshold)
                         obj.beads(i) = 0; 
                         break;
                     end
                 end
                     
          end
          
            %have to check that these guys are beads actually 
           num_beads = length(find(obj.beads==1))
           num_non_beads = length(find(obj.beads==0)) 
            
        end
        
        function obj = reorder_events_by_channel(obj, channel_name)
            
           channel =  obj.name_channel_map(channel_name);
           
           channel_data = obj.data(:,channel);
           [ ~, new_order] = sort(channel_data);
           obj.data = obj.data(new_order,:);
            
        end
        
        function [cdata_fracs] = split_object_into_fractions(obj, channel_name, num_fracs)
            
            cdata = obj.reorder_events_by_channel(channel_name);
            
            num_cells = floor(size(obj.data,1)/num_fracs);
            cdata_fracs = cell(1,num_fracs);
            
            for i=1:num_fracs
                
                cdata_fracs{i} = cdata;
                start_index = ((i-1)*num_cells)+1;
                finish_index = (i)*num_cells;
                cdata_fracs{i}.data = cdata.data(start_index:finish_index,:);
            end
            
        end
            
        

        



        
        function plot_events_vs_time(obj, channel_name, window_size, varargin)
            
           channel =  obj.name_channel_map(channel_name);
           optargin = size(varargin,2);
           num_events = size(obj.data,1);
           
           
           
           if(optargin > 0)
               
               time_data = varargin{1};
           else
               time_data = 1:num_events;
           end
           
           [ wanderlust_value , wanderlust_order] = sort(time_data);
            time_data = time_data(wanderlust_order);
            event_data = obj.data(wanderlust_order,channel);
           
            
           smooth_data = zeros(1,num_events-window_size);
           smooth_index = zeros(1,num_events-window_size);
            
           for i=1:num_events-window_size
            
               smooth_data(i) = median(event_data(i:i+window_size));
               %smooth_index(i) = median(event_data(i:i+window_size,1));
               %ordering_time = varargin{1};
               smooth_index(i) = median(time_data(i:i+window_size));
           end
            
           %h=plot(smooth_index, smooth_data, 'LineWidth',2.0);
           %set(h,'Color',varargin{2});
           
           plot(smooth_index, smooth_data,'LineWidth',2.0);
           
           %set(gca,'XTick',[]);
           xlabel('Wanderlust Order','FontSize',14);
           ystring = sprintf('%s level',channel_name);
           ylabel(ystring, 'FontSize',14);
            %scatter(event_data(:,1), event_data(:,channel));
            
        end
        
        function plot_bead_events(obj, channel, window_size)
            
           events = find(obj.beads ==1);
           obj.plot_events_vs_time(channel, window_size, events);
            
        end

        
         
         
        function [XT, YT] = interpolate_events(obj, total_events, channel_name, mode)
            
           channel = obj.name_channel_map(channel_name);
           [num_events,~] = size(obj.data);
           
           XI = zeros(total_events,1);
            
           for i=1:total_events
               
               XI(i) = i;
               
           end
           
           indices = obj.data(:,obj.eventnum_channel);
           XI(indices)=[];
            
           YI = interp1(obj.data(:,obj.eventnum_channel),obj.data(:,channel),XI, mode); 
           
           XT = zeros(total_events, 1);
           YT = zeros(total_events, 1);
           
           for i=1:num_events
               
               xval = obj.data(i,obj.eventnum_channel);
               yval = obj.data(i, channel);
               
               XT(xval) = xval;
               YT(xval) = yval;
              
           end
           
           for i=1:length(XI)
               
              
               
               XT(XI(i)) = XI(i);
               YT(XI(i)) = YI(i);
              
           end
           
           
           
        end
        

        function obj = compute_acquisition_rates(obj, window_size)
            
            num_events = length(obj.acquisition_times);
            obj.acquisition_rate = zeros(1, num_events-window_size);
            obj.acquisition_rate_index = zeros(1, num_events-window_size);
            half_window = window_size/2;
            
             
            for i=1:num_events-window_size
                start_time = obj.acquisition_times(i);
                end_time = obj.acquisition_times(i+window_size);
                duration = end_time - start_time;
                if(duration == 0)
                    duration = .00000001;
                end
                rate = window_size/duration;
                %number of acquisitions per millisecond
                
                obj.acquisition_rate(i) = rate;
                obj.acquisition_rate_index(i) = obj.acquisition_times(i+half_window);
            end  
            
        end
        
        function plot_acquisition_rate(obj, varargin)
            
            optargin = size(varargin, 2);
            
            if(optargin == 0)
                plot(obj.acquisition_rate_index, obj.acquisition_rate , 'r');
            else          
                plot(obj.acquisition_rate_index, obj.acquisition_rate , varargin{1});
            end
            
        end
        
        function plot_smooth_events(obj, channel_name, varargin)
            
            optargin = size(varargin, 2);
            channel = obj.name_channel_map(channel_name);
            if(optargin == 0)
                plot(obj.smooth_index, obj.smooth_data(:,channel));
            else          
                plot(obj.smooth_index, obj.smooth_data(:, channel), varargin{1});
            end
            
        end
        
        function obj = compute_smooth_corrs(obj)
            
            num_channels = size(obj.smooth_data,2);
            obj.smooth_corrs = zeros(num_channels, num_channels); 
            
            for i=1:num_channels
                for j=1:num_channels
                    obj.smooth_corrs(i,j) = corr(obj.smooth_data(:,i), obj.smooth_data(:,j));
                end
            end
            
        end
        
        function obj = compute_corrections(obj)
            
           num_data_points = size(obj.smooth_data,1); 
           num_channels = size(obj.smooth_data,2);
           obj.smooth_diffs = zeros(num_data_points-1, num_channels);
           obj.smooth_data_recon = zeros(num_data_points, num_channels);
           
           obj.smooth_fluctuation_vector = zeros(num_data_points-1, 1);
           
            for i=1:num_channels
              
                    obj.smooth_diffs(:,i) = diff(obj.smooth_data(:,i));
                    obj.smooth_fluctuation_vector = obj.smooth_fluctuation_vector + obj.smooth_diffs(:,i);
            end
            
            obj.smooth_fluctuation_vector = obj.smooth_fluctuation_vector .* (1/num_channels);
            %average the difference trends 
            
            for i=1:num_channels
                
                obj.smooth_data_recon(1,i) = obj.smooth_data(1,i);
                
                for j = 1:(num_data_points-1)
                    obj.smooth_data_recon(j+1,i) = obj.smooth_data_recon(j,i)+obj.smooth_fluctuation_vector(j);
                end
                
            end
            
        end
        
        function [correction_function] = create_correction(obj, base_level, acquisition_times)
           
           num_data_points = length(obj.acquisition_times); 
           correction_function = zeros(num_data_points, 1); 
           
           smooth_data_points = length(obj.smooth_index); 
            
           smooth_data_recon = zeros(smooth_data_points, 1);
           smooth_data_recon(1) = base_level;
           
           for i = 1:(smooth_data_points-1)
           
               smooth_data_recon(i+1) = smooth_data_recon(i) + obj.smooth_fluctuation_vector(i);
           
           end
           
           
           
           first_zero = -1;
           last_zero = 0;
           
           for i=2:length(acquisition_times)
            
                if(acquisition_times(i) == acquisition_times(i-1))   
                
                  
                    if(first_zero == -1)
                        first_zero = i-1;
                       
                    end
                
                else
                
                     if(first_zero > -1)
                  
                        last_zero = i-1;
                   
                        increments = (last_zero - first_zero)+1;
                        if(increments==0)
                            
                           sprintf('WARNING! Increments is 0 \n'); 
                        end
                        step = 1/increments;
                        sprintf('fixing acquisition time\n');
                        
                        for j=(first_zero+1):last_zero 
                      
                            acquisition_times(j) = acquisition_times(j) + step;
                            step = step+(1/increments);
                       
                        end
               
                        first_zero = -1;
                        last_zero = 0;
                   
                     end
               
                end
                   
           end    
           
          
           
           first_zero = -1;
           last_zero = 0;
           smooth_index = obj.smooth_index;
           
           for i=2:length(smooth_index)
            
                if(smooth_index(i) == smooth_index(i-1))   
                
                  
                    if(first_zero == -1)
                        first_zero = i-1;
                       
                    end
                
                else
                
                     if(first_zero > -1)
                  
                        last_zero = i-1;
                   
                        increments = (last_zero - first_zero)+1;
                        step = 1/increments;
                        sprintf('fixing acquisition time\n');
                        
                        for j=(first_zero+1):last_zero 
                      
                            smooth_index(j) = smooth_index(j) + step;
                            step = step+(1/increments);
                       
                        end
               
                        first_zero = -1;
                        last_zero = 0;
                   
                     end
               
                end
                   
           end  
           
           %if the last two timestamps are the same the above loop doesn't
           %fix them, so adding this here -rlf 20110811
           if(first_zero > -1)
               
               last_zero = i;
               
               increments = (last_zero - first_zero)+1;
               step = 1/increments;
               sprintf('fixing acquisition time\n');
               
               for j=(first_zero+1):last_zero
                   
                   smooth_index(j) = smooth_index(j) + step;
                   step = step+(1/increments);
                   
               end
               
           end
           
           
           full_data_recon = interp1(smooth_index, smooth_data_recon, acquisition_times, 'linear', 'extrap');
           
           
           correction_function = base_level - full_data_recon;
           
           
          
           
         end
        
          function obj = correct_channels(obj, sensitivity_threshold)
            %trying to return a modified version of this object
            
            num_channels = size(obj.data,2)
            %correcting all channels besides time and cell length
            
            
            
            
            for i = 3:num_channels
               
               channel = i;
               
               %obj.channel_averages = obj.channel_averages .* (1/length(obj.acquisition_times));
               
               base_level = median(obj.data(:,channel));
               
               correction_function = obj.create_correction(base_level, obj.data(:,1));
               
               low_sensitivity_indices = find(obj.data(:,channel)<sensitivity_threshold);
               
               correction_function(low_sensitivity_indices)=0;
               
               obj.data(:,channel) = obj.data(:,channel) + correction_function;
               
            end
            
%             num_channels = length(obj.dna_channels)
%             %marker channels, dna channels also 
%           
%             
%             for i = 1:num_channels
%                
%                channel = obj.dna_channels(i); 
%           
%                %obj.channel_averages = obj.channel_averages .* (1/length(obj.acquisition_times));
%                
%                base_level = median(obj.data(:,channel));
%                
%                correction_function = obj.create_correction(base_level, obj.data(:,1));
%                
%                obj.data(:,channel) = obj.data(:,channel) + correction_function;
%                
%             end
            
            
           
          end
        
          function obj = find_bead_data(obj)
             
             %obj = obj.identify_beads();
             %can add that back for automation purposes
             
             obj.acquisition_times = obj.data(:,1);
              
             num_beads = length(find(obj.beads==1))
             num_non_beads = length(find(obj.beads==0))
             
             bead_events = find(obj.beads==1);
             all_channels = [1 obj.bead_channels];
             obj.bead_data = obj.data(bead_events,all_channels);
                
              
          end
          
          function obj = bead_normalize(obj, corr_threshold, approx_percentage, sensitivity_threshold)
              
             
              obj = obj.identify_beads();
              sprintf('done identifying beads \n');
              
              obj = obj.find_bead_data();
              
              sprintf('seperated bead data \n');
              
              obj = obj.compute_smooth_events_window_opt(corr_threshold);
              
              sprintf('seperated smooth events \n');
              
              obj = obj.compute_corrections();
              
              sprintf('computed corrections \n');
              
              obj = obj.correct_channels(sensitivity_threshold);
              
              sprintf('corrected channels \n');
              
          end
       
          function obj = bead_normalize_fix_window(obj, window_size, sensitivity_threshold)
              
             
              %obj = obj.identify_beads();
              %sprintf('done identifying beads \n');
              
             
              
              obj = obj.find_bead_data();
              
              sprintf('seperated bead data \n')
              
              obj = obj.compute_smooth_events_vs_time(window_size);
              
              sprintf('seperated smooth events \n')
              
              obj = obj.compute_corrections();
              
              sprintf('computed corrections \n')
              obj = obj.correct_channels(sensitivity_threshold);
              
              sprintf('corrected channels \n')
              %obj = obj.remove_beads();
          end


       
        
    end
    
end
    


