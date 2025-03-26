classdef cytof_series
   
    properties  
        
        datasets
        datalabels
        series_name
        control
        marker_channels
    end
    
    methods
        
            
        
    function obj = cytof_series(datasets, datalabels, series_name)
    
        
        obj.datasets = datasets;
        obj.datalabels = datalabels;
        obj.series_name = series_name;
        cdata = datasets{1};
        obj.marker_channels = cdata.channel_name_map(cdata.marker_channels);
    end
    
    %%%%%%
    function obj = add_dataset(obj, cdata, datalabel)
        %datasets = obj.datasets;
        datalength = length(obj.datasets);
        obj.datasets{datalength+1} = cdata;
        obj.datalabels{datalength+1} = datalabel;
        
    end
    %%%%%%
    
    function obj = add_datasets(obj, cdatasets, datalabels)
       
        for i=1:length(cdatasets)
           
            obj = obj.add_dataset(cdatasets{i},datalabels{i});
            
        end
        
    end
    
    function cdataobj = make_combined_cytof_data_object(obj)
        
       cdataobj = obj.datasets{1};
       for i=2:length(obj.datasets{1})
           
           cdataobj = cdataobj.append_cdata_object(obj.dataests{i});
           
       end
        
    end
    
    
    
    function [mi_vector ] = compute_series_R(obj, channel1_name, channel2_name, timepoints)
        
        
       
        num_timepoints = length(timepoints); 
        mi_vector = zeros(1, num_timepoints);
      
     
      
     
         for i = 1:num_timepoints
                
              cdata = obj.datasets{timepoints(i)}; 
              R = cdata.corrcoef_edge(channel1_name,channel2_name);
              mi_vector(i) = R(2,1);
               

         end
         
         bar(mi_vector,'k');
         %bar(pval_vector,'g');
         set(gca,'XTick',[]);
         set(gca,'YTick',[]);
         set(gca, 'XTickLabel','');
         set(gca, 'YTickLabel','');
         
    end
    
    
    function [mi_vector, pval_vector] = compute_series_mi(obj, channel1_name, channel2_name, timepoints, threshold, varargin)
        
        
       
        
        num_timepoints = length(timepoints); 
        mi_vector = zeros(1, num_timepoints);
        pval_vector = zeros(1, num_timepoints);
        maxy = 0;
        
        
        
        
        for i=1:num_timepoints
             
              cdata = obj.datasets{timepoints(i)};
              
              
              [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdata, channel1_name, channel2_name, 50, 255);
              
              if(maxy1>maxy) 
                  
                    maxy = maxy1;
              end
            
            
        end
            
       for i=1:length(varargin)
           
          if(strcmp(varargin{i},'maxy'))
              
              maxy = varargin{i+1};
          end
          
       end
        
        
        
         for i = 1:num_timepoints
                
               
                [mi_vector(i), pval_vector(i)] = obj.compute_mi(channel1_name,channel2_name,timepoints(i), 'maxy', maxy, 'noise_threshold', threshold);
                %[mi_vector(i), pval_vector(i)] = obj.compute_mi(channel1_name,channel2_name,timepoints(i), 'noise_threshold', .8);

         end
         
%          bar(mi_vector,'k');
%          %bar(pval_vector,'g');
%          set(gca,'XTick',[]);
%          set(gca,'YTick',[]);
%          set(gca, 'XTickLabel','');
%          set(gca, 'YTickLabel','');
         
    end
    

  
   
    
   function obj = gate_high_low_series(obj, channel_name, percent_start, percent_end)
        
       for i=1:length(obj.datasets)
          
           cdata = obj.datasets{i};
           
           cdata = cdata.high_low_gate_events(channel_name, percent_start, percent_end);
           
           obj.datasets{i} = cdata;
           
       end
        
        
        
   end
   
   function obj = manual_gate_series(obj, channel1_name, channel2_name)
       
            fhandle = figure;
           
            cdata = obj.datasets{1};
            
            channel1 = cdata.name_channel_map(channel1_name);
            channel2 = cdata.name_channel_map(channel2_name);
            
            channel_data = [cdata.data(:,channel1) cdata.data(:,channel2)]; 
             [~, density, x, y] = kde2d([channel_data(:,1) channel_data(:,2)], 256);

            plot(channel_data(:,1),channel_data(:,2),'b.','MarkerSize',5)
            hold on,
            contour(x, y, density, 30);
            
            
           
           
        
     
        
            rect = getrect(fhandle)
            left = rect(1)
            bottom = rect(2)
            width = rect(3)
            height = rect(4)
       
        
            for i=1:length(obj.datasets)
                
               cdata = obj.datasets{i};
               cdata = cdata.manual_gate_events_box(channel1_name, channel2_name, left, bottom, width, height);
               obj.datasets{i} = cdata;
                
                
            end
        
        
   end 
   
    function [csets] = conditional_gated_series(obj, channel_name, num_conditions)
        
       
        maxval = 0;
        minval = 0;
       
       
       conditional_series = cell(1,num_conditions);
       conditional_series_label = cell(1,num_conditions);
       condition_increment = 1/num_conditions;
       
       for j = 1:num_conditions
           
           new_datasets = cell(1, length(obj.datasets));
           
           for i=1:length(obj.datasets)
        
                
               percent_start = 0+ (condition_increment*(j-1));
               percent_end = percent_start+condition_increment;
               cdata = obj.datasets{i};
               cdata = cdata.high_low_gate_events(channel_name, percent_start, percent_end);
               new_datasets{i} = cdata;
              
               
             
               
           end
           
           s = sprintf('%s_conditional_series_%d',channel_name, j);
           cseries = cytof_series(new_datasets, obj.datalabels, s);
          
           conditional_series{j} = cseries; 
           
           conditional_series_label{j} = s;
       end
   
        
      
       csets = cytof_sets(conditional_series, conditional_series_label, 'conditional_set');
        
        
    end
    
    function [cseries] = conditional_interval_gated_series(obj, channel_name, intervals, timepoint)
        
    
      
       
       num_conditions = length(intervals)-1; 
       new_datasets = cell(1, num_conditions);
       conditional_series_label = cell(1,num_conditions);
       
       for j = 1:num_conditions
           
           
           cdata = obj.datasets{timepoint};
          
           
           cdata = cdata.high_low_gate_events_interval(channel_name, intervals(j), intervals(j+1))
           new_datasets{j} = cdata;
           conditional_series_label{j} = sprintf('interval_%d',j);
           
               
       end
           
       s = sprintf('%s_conditional_series_%d',channel_name, j);
       cseries = cytof_series(new_datasets, conditional_series_label, s);
          
        
    end
    
    function [ data1_bar ] = level_bargraph(obj, channel1_name)
        
        num_datasets = length(obj.datasets)
        data1_bar = zeros(num_datasets,1);
      
        for j=1:num_datasets
            
            cdataj = obj.datasets{j};
          
            data1_bar(j) = cdataj.get_channel_mean(channel1_name);
           
        end
        
        bar(data1_bar);
        
    end
    
     function [bar_data] = median_level_image(obj, marker_channels,timepoints)
            
           
            num_channels = length(marker_channels);
            bar_data = [];
           
            channel_names = cell(1,num_channels);
            for j=1:num_channels
              cdata = obj.datasets{1};  
              %channel_name = cdata.channel_name_map{marker_channels(j)}
              channel_name = marker_channels{j};
              series_data = transpose(obj.level_bargraph(channel_name));
%                for k=1:length(series_data)
%                   
%                    series_data(k) = (series_data(k)-series_data(1))/series_data(1)
%                     if(series_data(k)<1)
%                         series_data(k) = -1/series_data(k);
% %                    end
%                end
              
              bar_data = [bar_data; series_data(timepoints)];
              channel_names{j} = channel_name;
            end
%             max_val1 = max(max(bar_data));
%             max_val2 = abs(min(min(bar_data)));
%             max_val = max(max_val1, max_val2);
%             
%             clim = [-max_val, max_val];
             load continuous_BuPu9.mat
             colormap(continuous_BuPu9)
             imagesc(bar_data,[0, 6]);
             channel_names
              set(gca,'ytick',1:length(channel_names));
             %set(gca,'YTick',[]);
             set(gca,'yticklabel',channel_names);
             set(gca,'XTick',[]);
             %axis off
          %  bar(bar_data);
      end
            
    
    
    function [data1_bar data2_bar] = pairwise_bargraph_datasets(obj, cseries2, channel1_name)
        
        
        num_datasets = length(obj.datasets)
        data1_bar = zeros(num_datasets,1);
        data2_bar = zeros(num_datasets,1);
        
        
        for j=1:num_datasets
            
            cdataj = obj.datasets{j};
            cdatak = cseries2.datasets{j};
            data1_bar(j) = cdataj.get_channel_mean(channel1_name);
            data2_bar(j) = cdatak.get_channel_mean(channel1_name);
        
        end
    
        bar([data1_bar data2_bar]);
%         datax = linspace(1,num_datasets,num_datasets);
%         dataxi = linspace(1,num_datasets,200);
%         datay1=interp1(datax,data1_bar,dataxi,'cubic');
%         datay2=interp1(datax,data2_bar,dataxi,'cubic');
%         plot(dataxi,datay1,'r','LineWidth',3.0);
%         hold
%         plot(dataxi,datay2,'g','LineWidth',3.0);
        
       
        
    end
       
    
    function [data1_bar data2_bar] = edge_activity_bargraph(obj, cseries2, channel1_name, channel2_name)
        
        
        num_datasets = length(obj.datasets)
        data1_bar = zeros(num_datasets,1);
        data2_bar = zeros(num_datasets,1);
        
        
        for j=1:num_datasets
            
            
            
            data1_bar(j) = obj.compute_activity(channel1_name,channel2_name,j);
            %data2_bar(j) = cseries2.compute_activity(channel1_name,channel2_name,j);
        
        end
    
        %bar([data1_bar data2_bar]);
        bar(data1_bar,'c');
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
    end
    
    function [mi1_bar mi2_bar pval1_bar pval2_bar] = edge_mi_bargraph(obj, cseries2, channel1_name, channel2_name)
        
        
        num_datasets = length(obj.datasets)
        mi1_bar = zeros(num_datasets,1);
        mi2_bar = zeros(num_datasets,1);
        
        
        for j=1:num_datasets
            
            
            
            [mi1_bar(j), pval1_bar(j)] = obj.compute_mi(channel1_name,channel2_name,j);
            [mi2_bar(j), pval2_bar(j)] = cseries2.compute_mi(channel1_name,channel2_name,j);
        
        end
    
        bar([mi1_bar mi2_bar]);
        
    end
    
    
    function [ce1_bar ce2_bar] = edge_centropy_bargraph(obj, cseries2, channel1_name, channel2_name)
        
        
        num_datasets = min(length(cseries2.datasets),length(obj.datasets));
        ce1_bar = zeros(num_datasets,1);
        ce2_bar = zeros(num_datasets,1);
        
        
        for j=1:num_datasets
            
            [ce1_bar(j),ce2_bar(j)] = obj.compute_centropy_compare(cseries2, channel1_name, channel2_name, j)
          
        
        end
    
        bar([ce1_bar ce2_bar]);
        
    end
    
    
   
   
    
    function pairwise_visualize_compare(obj,cseries2,channel1_name, channel2_name, time_point)
       
        cdataj = obj.datasets{time_point};
        cdatak = cseries2.datasets{time_point};
        
         [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
         [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
                minx = max(minx1,minx2);
                miny = max(miny1, miny2);
                maxx = min(maxx1, maxx2);
                maxy = min(maxy1, maxy2);
                
         subplot(1,2,1);       
            
         [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name,'Limits', [minx miny maxx maxy]);

        
         subplot(1,2,2);
               
         [points_x, points_y] = cdatak.pairwise_visualize(channel1_name,channel2_name,'Limits', [minx miny maxx maxy]);
        
        
        
    end
    
    function pairwise_correlation_visualization_datasets_compare(obj, cseries2, channel1_name, channel2_name, time_points,varargin)
    
       
       num_plot_rows = 2;
       fix_limits = 0;
       set_maxy = 0;
       maxy_val = 0;
       
       %if(show_density) 
       %    num_plot_rows = num_plot_rows*2;
       %end
       
   
       num_plot_cols = length(time_points);
       
       for i=1:length(varargin)
           
          if(strcmp(varargin{i},'fix_limits'))
              
             fix_limits = 1; 
          end
          
          if(strcmp(varargin{i},'maxy'))
              
             set_maxy = 1;
             maxy_val = varargin
          end
           
         
       end
       
       
       
    
       
       maxy = 0;
     
       for j = 1:length(time_points)
          
            time_point = time_points(j);

            subplot(num_plot_rows, num_plot_cols, j);
            cdataj = obj.datasets{time_point};
           
           
           
            channel1 = cdataj.name_channel_map(channel1_name);
            channel2 = cdataj.name_channel_map(channel2_name);
            cdatak = cseries2.datasets{time_point};
            [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
            [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(cdatak, channel1_name, channel2_name, 50, 255);
            
            if(maxy1>maxy) 
                maxy = maxy1;
            end
            if(maxy2>maxy)
                maxy = maxy2;
            end
            
          
       end
       
        
       
       for j = 1:length(time_points)
          
            time_point = time_points(j);

            subplot(num_plot_rows, num_plot_cols, j);
            cdataj = obj.datasets{time_point};
           
           
           
            channel1 = cdataj.name_channel_map(channel1_name);
            channel2 = cdataj.name_channel_map(channel2_name);
            cdatak = cseries2.datasets{time_point};
            
            if(fix_limits ==0)
                
               [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name);

        
               subplot(num_plot_rows, num_plot_cols, j+num_plot_cols);
               
               [points_x, points_y] = cdatak.pairwise_visualize(channel1_name,channel2_name);

            else
              
                [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
                [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(cdatak, channel1_name, channel2_name, 50, 255);
                
                maxx = min(maxx1, maxx2);
                
            
               [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name,'Limits',[0 0 maxx maxy]);
               %xlim([0 maxx]);
               %ylim([0 maxy]);
               set(gca,'XTick',[1:5]);
               set(gca,'YTick',[1:12])
        
               subplot(num_plot_rows, num_plot_cols, j+num_plot_cols);
               
               [points_x, points_y] = cdatak.pairwise_visualize(channel1_name,channel2_name,'Limits',[0 0 maxx maxy]);
                
               %xlim([0 maxx]);
               %ylim([0 maxy]);
               set(gca,'XTick',[1:5]);
               set(gca,'YTick',[1:12])
            end
       end
       
         
       
    end
    
    function write_series_edge_csv(obj, channel1_name, channel2_name, filename, timepoints)
        
        
        num_timepoints = length(timepoints);
        cd = obj.datasets{timepoints(1)};
        min_points = size(cd.data,1);
        for i=1:num_timepoints
            cdata = obj.datasets{timepoints(i)};
            
            if(min_points>size(cdata.data,1))
                min_points = size(cdata.data,1);
            end
            
        end
        
        %min_points = 7845;
        %m = zeros(num_timepoints*2,max_points);
        
        for i=1:num_timepoints
            cdata = obj.datasets{timepoints(i)};
            channel1 = cdata.name_channel_map(channel1_name);
            channel2 = cdata.name_channel_map(channel2_name);
            
           
            num_datapoints = size(cdata.data,1);
            index1 = (i*2)-1;
            index2 = (i*2);
            
            %m = [transpose(cdata.data(:,channel1));transpose(cdata.data(:,channel2))];
            m = [transpose(cdata.data(1:min_points,channel1));transpose(cdata.data(1:min_points,channel2))];
                
            %m(index1,1:num_datapoints) = transpose(cdata.data(:,channel1));
            %m(index2,1:num_datapoints) = transpose(cdata.data(:,channel2));
            %filename = sprintf('%s_t%d.csv',filename_root,i);
            %csvwrite(filename,m);
            dlmwrite(filename,m,'-append','delimiter',',');
        end
        
        %csvwrite(filename,m);
        
    end
    
    function pairwise_correlation_visualization_datasets(obj, channel1_name, channel2_name, time_points,varargin)
    
       
       num_plot_rows = 1;
       fix_limits = 0;
       scatter = 0;
       kdeplot = 0;
       regression = 0;
       axes_specified = 0;
       axes_handle = 0;
       %if(show_density) 
       %    num_plot_rows = num_plot_rows*2;
       %end
       
        set_maxy = 0;
        maxy_val = 0;
       
       num_plot_cols = length(time_points);
       
       for i=1:length(varargin)
           
           if(strcmp(varargin{i},'axes'))
               
              axes_handle = varargin{i+1}; 
           end
           
          if(strcmp(varargin{i},'fix_limits'))
              
             fix_limits = 1; 
          end
          if(strcmp(varargin{i},'scatter'))
              
             scatter = 1;
             
          end
           
          if(strcmp(varargin{i},'kde'))
              
             kdeplot = 1;
             
          end
          
          if(strcmp(varargin{i},'sigmoid'))
              
             regression = 1;
             
          end

          if(strcmp(varargin{i},'linear'))
              
             regression = 2;
             
          end
          
          if(strcmp(varargin{i},'condmean'))
              
             regression = 3;
             
          end
          
         
          
           if(strcmp(varargin{i},'maxy'))
              
             set_maxy = 1;
             maxy_val = varargin{i+1};
          end
          
       end
       
       
       
    
       
       maxy = 0;
       minx = 100;
     
       for j = 1:length(time_points)
          
            time_point = time_points(j);

          
            cdataj = obj.datasets{time_point};
           
           
           
            channel1 = cdataj.name_channel_map(channel1_name);
            channel2 = cdataj.name_channel_map(channel2_name);
           
            [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
         
            if(maxy1>maxy) 
                maxy = maxy1;
            end
           
            if(minx1<minx)
                minx = minx1;
            end    
       end
       
       %overwrite it if its fixed.
       if(set_maxy==1)
           maxy = maxy_val;
       end
       
       for j = 1:length(time_points)
          
            time_point = time_points(j);

            subplot(num_plot_rows, num_plot_cols, j);
           %subplot(num_plot_cols, num_plot_rows, j);
            cdataj = obj.datasets{time_point};
           
           
           
            channel1 = cdataj.name_channel_map(channel1_name);
            channel2 = cdataj.name_channel_map(channel2_name);
           
            
            if(fix_limits ==0)
                
               
               
               if(scatter == 1)
                   
                   cdataj.plot_2d_channel_scatter(channel1_name, channel2_name);
               elseif(kdeplot ==1 )
                   
                    cdataj.plot_2d_channel_density(channel1_name, channel2_name,'imagesc');
               else
                   [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name);
               end

            else
              
               if(scatter == 1)
                   
                   cdataj.plot_2d_channel_scatter(channel1_name, channel2_name);
               elseif(kdeplot ==1 )
                   
                    cdataj.plot_2d_channel_density(channel1_name, channel2_name,'imagesc');
               else
                   
                    
                    
                    [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdataj, channel1_name, channel2_name, 50, 255);
                    min_use_x = max([0, minx1]);
                    min_use_y = max([0, miny1]);
                    %limits = [0 0 maxx1 maxy];
                    
                    if(regression ==0)
                        
                      if(axes_specified ==0)     
                           [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name,'Limits',[min_use_x min_use_y maxx1 maxy]);
                      else
                           [points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name,'Limits',[min_use_x min_use_y maxx1 maxy],'axes',axes_handle);
                      end
                    end
                    
                    
                    if(regression ==1 )
                        [points_x, points_y] = cdataj.sigmoid_fit_edge(channel1_name,channel2_name,'Limits',[min_use_x min_use_y maxx1 maxy]);
                    end
                    
                    if(regression ==2)
                         [points_x, points_y] = cdataj.linear_fit_edge(channel1_name,channel2_name,'Limits',[min_use_x min_use_y maxx1 maxy] );
                    end
                    
                    if(regression == 3)
                        %disp('doing');
                        [x, y] = cdataj.conditional_mean_edge(channel1_name, channel2_name,'Limits',[min_use_x min_use_y maxx1 maxy] )
                        
                    end
                    
                    
                   %[points_x, points_y] = cdataj.pairwise_visualize(channel1_name,channel2_name,'MinMaxFitX','MinMaxY',0, maxy);
%                     if(j==1)
%                      
%                     else
%                     [points_x, points_y] = cdataj.sigmoid_fit_edge(channel1_name,channel2_name);
%                     end
%                     [responding_fractions] = cdataj.pairwise_visualize_responding_fraction(channel1_name,channel2_name, 8, 0.6806,'Limits',[0 0 maxx1 maxy]);
                    %cdataj.plot_2d_channel_density(channel1_name, channel2_name, 'limits',[0 0 maxx1 maxy], 'imagesc');
%                     cdataj.plot_2d_channel_density(channel1_name, channel2_name, 'imagesc');
%                     line_xs = linspace(1,7,7);
%                     xlim = get(gca,'XLim')
%                     line_xs = line_xs.*((xlim(2)-xlim(1))/8)
%                     plot_as_vertical_lines(line_xs,'w');
                    %cdataj.plot_2d_channel_scatter(channel1_name, channel2_name);
               end

                
                
                
                    channel1_name = upper(channel1_name)
                    channel2_name = upper(channel2_name)
                    channel1_name = strrep(extractAfter(channel1_name, '_'), '_', '');
                    channel2_name = strrep(extractAfter(channel2_name, '_'), '_', '');
                    xlabel(channel1_name,'FontSize',16);
                    ylabel(channel2_name,'FontSize',16);
                    set(gca,'TickLabelInterpreter','none');
                    title("Allah", 'FontSize',16);
        
               
            end
       end
       
         
       
    end
    
    
    
    function [delta_entropy1, delta_entropy2] = time_series_delta_compare(obj,cseries2, channel1_name, channel2_name, num_slicesx, num_slicesy)
        
        [delta_entropy1] =  obj.time_series_delta(channel1_name, channel2_name, num_slicesx, num_slicesy);
        [delta_entropy2] = cseries2.time_series_delta(channel1_name, channel2_name, num_slicesx, num_slicesy);
        bar([delta_entropy1 delta_entropy2]);
        
        
    
    
    end
    
    
    function [mi_diff_matrix] = pairwise_mi_timepoint_compare(obj, cseries2, channel_names, num_slicesx, num_slicesy,timepoint)
        
        mi1_matrix = obj.pairwise_mi_timepoint(channel_names, num_slicesx, num_slicesy, timepoint, 'no_plot');
        mi2_matrix = cseries2.pairwise_mi_timepoint(channel_names, num_slicesx, num_slicesy, timepoint, 'no_plot');
        mi_diff_matrix = mi1_matrix-mi2_matrix;
        CLIM = [0 0.04]; 
        colormap(redgreencmap); 
        imagesc(mi_diff_matrix, CLIM);

        set(gca,'ytick',1:length(channel_names));
        set(gca,'yticklabel',channel_names);
        xticklabel_rotate([1:length(channel_names)],45,channel_names);
         
        
    end
    
    
    
   
    
    function [mi, pval] = compute_mi(obj, channel1_name, channel2_name, timepoint, varargin)
                
               cdata = obj.datasets{timepoint};
        
               
                num_permutations = 1000;
                maxy = 0; 
                noise_threshold = 0.8;
                
                for i=1:length(varargin)
                 
                 
                   if(strcmp(varargin{i}, 'maxy'))
                    maxy = varargin{i+1};
                   end
                  
                   if(strcmp(varargin{i}, 'noise_threshold'))
                       noise_threshold = varargin{i+1};
                   end
                       
                      
                  
                end
      
                
               pval = 0; 
                
                
                if(maxy>0)     
                
                    %[mi, pval] = cdata.compute_dremi(channel1_name, channel2_name, .8, 'compute_pvalue', 10000);
                    [mi, pval] = cdata.compute_dremi(channel1_name, channel2_name, noise_threshold, 'maxy', maxy);
                else
                    
                    %[mi, pval] = cdata.compute_dremi(channel1_name, channel2_name, .8, 'compute_pvalue', 10000, 'maxy', maxy); 
                    [mi, pval] = cdata.compute_dremi(channel1_name, channel2_name, noise_threshold); 
                end
        
        
        
    end
    
    
     function [c] = compute_centropy(obj, channel1_name, channel2_name, timepoint)
                
               cdata = obj.datasets{timepoint};
        
                channel1 = cdata.name_channel_map(channel1_name);
                channel1_data = cdata.data(:,channel1);
                channel2 = cdata.name_channel_map(channel2_name);
                channel2_data = cdata.data(:,channel2);
      
                num_slicesx = 8;
                num_slicesy = 32;
        
                [minx, miny, maxx, maxy] = find_data_cutoffs(cdata, channel1_name, channel2_name, 50, 255);

        
                [~, ~, c] = delta_entropyreweight(cdata, channel1_name, channel2_name, num_slicesx, num_slicesy, 0, 0, maxx, maxy);
                %[mi, pval] = compute_dremi_pvalue(cdata, channel1_name, channel2_name, num_slicesx, num_slicesy, 1000);

        
        
        
     end
    
    function [c1,c2] = compute_centropy_compare(obj, cseries2, channel1_name, channel2_name, timepoint)
                
               cdata = obj.datasets{timepoint};
               cdata2 = cseries2.datasets{timepoint};
                channel1 = cdata.name_channel_map(channel1_name);
                channel1_data = cdata.data(:,channel1);
                channel2 = cdata.name_channel_map(channel2_name);
                channel2_data = cdata.data(:,channel2);
      
                num_slicesx = 8;
                num_slicesy = 8;
        
                [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(cdata, channel1_name, channel2_name, 50, 255);
                [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(cdata2, channel1_name, channel2_name, 50, 255);
                maxx = max(maxx1,maxx2);
                maxy = max(maxy1, maxy2);
                
                min_use_y = min(miny1, miny2);
                min_use_x = min(minx1, minx2);
        
                [~, ~, c1] = delta_entropyreweight(cdata, channel1_name, channel2_name, num_slicesx, num_slicesy, min_use_x, min_use_y, maxx, maxy);
                 [~, ~, c2] = delta_entropyreweight(cdata2, channel1_name, channel2_name, num_slicesx, num_slicesy, min_use_x, min_use_y, maxx, maxy);
                %[mi, pval] = compute_dremi_pvalue(cdata, channel1_name, channel2_name, num_slicesx, num_slicesy, 1000);

        
        
        
    end
    
    
    
    
    
    
    
    function [delta_entropy] =  time_series_delta(obj, channel1_name, channel2_name, num_slicesx, num_slicesy,varargin)
   


        delta_entropy = [];
        entropyy = [];
        conditional_entropy = [];

      
    
    
            for i=1:length(obj.datasets)
    
                cdata = obj.datasets{i};
                channel1 = cdata.name_channel_map(channel1_name);
                channel1_data = cdata.data(:,channel1);
                channel2 = cdata.name_channel_map(channel2_name);
                channel2_data = cdata.data(:,channel2);
      
        
                [minx, miny, maxx, maxy] = find_data_cutoffs(cdata, channel1_name, channel2_name, 50, 255);

        
                [d, e, c] = delta_entropyreweight(cdata, channel1_name, channel2_name, num_slicesx, num_slicesy, miny, miny, maxx, maxy);
                %[points_x, points_y, normalized_density] = cdata.pairwise_visualize(channel1_name,channel2_name,'no_plot');
                %[d, e, c] = delta_entropyreweight_density(normalized_density);
                %[d, e, c] = delta_entropyreweight_peakpts(points_x, points_y, num_slicesx, num_slicesy)
   
                
                delta_entropy = [delta_entropy; d];
                entropyy = [entropyy; e];
                conditional_entropy = [conditional_entropy; c];
        
      
        
            end
    
            entropyy
            conditional_entropy
            delta_entropy
            plotdata = delta_entropy;
    
            %figure('visible', 'on');
            optargin = length(varargin);
            if(optargin>0)
                if(strcmp(varargin{1},'plot'))
                    
                    bar(plotdata);
                end
            end
    
    
    end
    
    function [total_differential,mi1, mi2] = time_series_mi_compare(obj, cseries2, channel1_name, channel2_name, num_slicesx, num_slicesy)
 
        num_time_points = length(obj.datasets);
        mi1 = zeros(1,num_time_points);
        mi2 = zeros(1,num_time_points);
    
        mi1 = obj.time_series_delta(channel1_name, channel2_name,num_slicesx, num_slicesy);
        mi2 = cseries2.time_series_delta(channel1_name, channel2_name,num_slicesx, num_slicesy);
        total_differential = norm(mi1-mi2,1);

    end
    
    function [total_differentials] = compute_total_differentials(obj, cseries2, channel1_names, num_slicesx, num_slicesy)
        
        
       total_differentials = zeros(length(channel1_names),length(channel1_names));

        for i=1:length(channel1_names)
            for j = 1:length(channel1_names)
                if(i==j)
                    continue;
                end
                %filename1 = sprintf('F_B6_%s_%s_centropy',channel1_names{j}, channel2_names{j});
                %filename2 = sprintf('F_nod_%s_%s_centropy',channel1_names{j}, channel2_names{j});
                [total_differentials(i,j),~, ~] = obj.time_series_mi_compare(cseries2, channel1_names{i}, channel1_names{j}, num_slicesx, num_slicesy);
        
            end
        end
        
        imagesc(total_differentials)

        set(gca,'ytick',1:length(channel1_names));
        set(gca,'yticklabel',channel1_names);
        xticklabel_rotate([1:length(channel1_names)],45,channel1_names);
        
    end
    
    
    
     function [sensitivity_matrix] = compute_activity_over_time(obj, channel1_name, channel_names, channel1_axis, varargin)
        
        num_timepoints = length(obj.datasets);
        num_channels = length(channel_names);
        sensitivity_matrix = zeros(num_channels, num_timepoints); 
        
        for i=1:num_timepoints
        
        
            
           
            for j=1:num_channels
                
               
                
                channel2_name = channel_names{j};
                if(strcmp(channel2_name, channel1_name))
                    continue;
                end
                
                if(strcmp(channel1_axis,'x'))
    
                    sensitivity_matrix(j,i) = obj.compute_activity(channel1_name, channel2_name, i);
                    
                else
                   
                    sensitivity_matrix(j,i) = obj.compute_activity(channel2_name, channel1_name, i);
                    
                end
                
                    
                
            end
        end
        
        minval = min(min(sensitivity_matrix));
        maxval = max(max(sensitivity_matrix));
        
        if(length(varargin)>0)
            minval = varargin{1};
        elseif(length(varargin)>1)
            maxval = varargin{2};
        end
       
        
        colormap(pink);
        CLIM = [0 maxval];
        imagesc(sensitivity_matrix, CLIM);
        set(gca,'ytick',1:length(channel_names));
        set(gca,'yticklabel',channel_names);
        xticklabel_rotate([1:length(obj.datalabels)],45,obj.datalabels);

        
        
        
    end
    
    
    
    function [sensitivity_diff_matrix] = compute_activity_difference_over_time(obj, cseries2, channel1_name, channel_names, channel1_axis, varargin)
        
        num_timepoints = length(obj.datasets);
        num_channels = length(channel_names);
        sensitivity_diff_matrix = zeros(num_channels, num_timepoints); 
        
        for i=1:num_timepoints
        
        
            
           
            for j=1:num_channels
                
               
                
                channel2_name = channel_names{j};
                if(strcmp(channel2_name, channel1_name))
                    continue;
                end
                
                if(strcmp(channel1_axis,'x'))
    
                    sensitivity_diff_matrix(j,i) = obj.compute_activity_difference(cseries2, channel1_name, channel2_name, i);
                    
                else
                   
                    sensitivity_diff_matrix(j,i) = obj.compute_activity_difference(cseries2, channel2_name, channel1_name, i);
                    
                end
                
                    
                
            end
        end
        
       minval = min(min(sensitivity_diff_matrix));
       maxval = max(max(sensitivity_diff_matrix));
        
        if(length(varargin)>0)
            minval = varargin{1};
        elseif(length(varargin)>1)
            maxval = varargin{2};
        end
       
        if(abs(minval)>abs(maxval))
            maxval = -1*minval;
        else
            minval = -1*maxval;
        end
        
        colormap(redgreencmap)
        CLIM = [minval maxval];
        imagesc(sensitivity_diff_matrix, CLIM);
        set(gca,'ytick',1:length(channel_names));
        set(gca,'yticklabel',channel_names);
        xticklabel_rotate([1:length(obj.datalabels)],45,obj.datalabels);

        
        
        
    end
    
    
    function [sensitivity, area1, area2] = compute_activity_difference(obj, cseries2, channel1_name, channel2_name, timepoint, varargin)
    
        cdata1 = obj.datasets{timepoint};
        cdata2 = cseries2.datasets{timepoint};
        [points_x1, points_y1] = cdata1.pairwise_visualize(channel1_name,channel2_name,'no_plot');
        [points_x2, points_y2] = cdata2.pairwise_visualize(channel1_name,channel2_name,'no_plot');

            
        range_min = max(min(points_x1),min(points_x2));
        range_max = min(max(points_x1),max(points_x2));
   
        estimation_pts = linspace(range_min, range_max, 100);
        [points_y1_estimate] = interp1(points_x1, points_y1, estimation_pts, 'linear');
        [points_y2_estimate] = interp1(points_x2, points_y2, estimation_pts, 'linear');
   
        area1 = norm(points_y1_estimate,1);
        area2 = norm(points_y2_estimate,1);
        
        delta = abs(estimation_pts(1)-estimation_pts(2))*10000;
        area1 = area1*delta;
        area2 = area2*delta;
        
        
        
%         if(area1<area2) 
%             sensitivity = -1;
%         else
%             sensitivity = 1;
%             
%         end
          sensitivity = (area1-area2)./area1;

         if(length(varargin)>0)
            num_args = length(varargin);
            
            for i=1:length(varargin)
        
             
               if(strcmp(varargin{i},'mi_threshold'))
                   
                  
                  [mi1,pval1] = obj.compute_mi(channel1_name, channel2_name, timepoint);
                  [mi2, pval2] = cseries2.compute_mi(channel1_name, channel2_name, timepoint);
                  if((pval1>0.05)&(pval2>0.05))
                      sensitivity = 0;
                  end
                  
                  
                   
               end
            end
            
           
         end
        
         
        
         
    end

    function [auc_values] = compute_series_auc(obj, channel1_name, channel2_name, timepoints)
        
        auc_values = zeros(1,length(timepoints));
       for i=1:length(timepoints)
           
          auc_values(i) = obj.compute_activity(channel1_name, channel2_name, timepoints(i));
           
       end
        
    end
    
    function [area1] = compute_activity(obj, channel1_name, channel2_name, timepoint, varargin)
    
        cdata1 = obj.datasets{timepoint};
      
        [points_x1, points_y1] = cdata1.pairwise_visualize(channel1_name,channel2_name,'no_plot');
       

            
        range_min = min(points_x1);
        range_max = max(points_x1);
   
        estimation_pts = linspace(range_min, range_max, 100);
        [points_y1_estimate] = interp1(points_x1, points_y1, estimation_pts, 'linear');
        
        delta = abs(estimation_pts(1)-estimation_pts(2));
   
        area1 = norm(points_y1_estimate,1);
        area1 = area1*delta;
            
        if(length(varargin)>0)
            num_args = length(varargin);
            
            for i=1:length(varargin)
        
               if(strcmp(varargin{i}, 'plot'))
                    plot(estimation_pts, points_y1_estimate,'r.');
                
               end
               
               
               if(strcmp(varargin{i},'mi_weight'))
                   
                  threshold = varargin{i+1}; 
                  mi = obj.compute_mi(channel1_name,channel2_name, timepoint);
                  if(mi<threshold)
                      area1 = area1*mi;
                  else
                      area1 = area1*1;
                  end
                  
                   
               end
               
               if(strcmp(varargin{i},'mi_threshold'))
                   
                  
                  [mi,pval] = obj.compute_mi(channel1_name,channel2_name, timepoint);
                  if(pval>0)
                      area1 = 0;
                  else
                      area1 = area1;
                  end
                  
                   
               end
               
               
               
            end
            
           
        end
    end
    
   
    function [sensitivity_diff_matrix] = compute_pairwise_activity_difference(obj, cseries2, channel1_names, varargin)
        
        num_channels = length(channel_names);
        sensitivity_diff_matrix = zeros(num_channels, num_channels);
        num_time_points = length(obj.datasets);
        
        for i=1:num_channels
            for j=1:num_channels
                
                if(i==j) 
                    continue;
                end
                
                channel1_name = channel_names{i};
                channel2_name = channel_names{j};
                for k=1:num_time_points
                    sensitivity_diff_matrix(i,j) = sensitivity_diff_matrix(i,j) + obj.compute_activity_difference(cseries2, channel1_name, channel2_name, k);
                end
                    
            end
        end
        minval = min(min(sensitivity_diff_matrix));
        maxval = max(max(sensitivity_diff_matrix));
        
        if(length(varargin)>0)
            minval = varargin{1};
        elseif(length(varargin)>1)
            maxval = varargin{2};
        end
        
        if(abs(minval)>abs(maxval))
            maxval = -1*minval;
        else
            minval = -1*maxval;
        end
        
        
        
        
        colormap(redgreencmap)
        CLIM = [minval, maxval];
        imagesc(sensitivity_diff_matrix, CLIM);
        set(gca,'ytick',1:length(channel_names));
        set(gca,'yticklabel',channel_names);
        xticklabel_rotate([1:length(channel_names)],45,channel_names);

        
    end
    
    
    function [activity_differential] = gen_all_timepoint_pairwise_activity_difference(obj, cseries2, channel_names, channel_labels)
        
         num_timepoints = length(obj.datasets);
         activity_differential = cell(1,num_timepoints);
         for i=1:num_timepoints
             
             activity_differential{i} = obj.compute_pairwise_activity_difference_timepoint(cseries2, channel_names, channel_labels, i); 
             
         end
        
        
        
    end
    
    function [activity_matrices1, activity_matrices2] = gen_all_timepoint_pairwise_activities(obj, cseries2, channel_names, channel_labels)
        
        num_timepoints = length(obj.datasets);
         activity_matrices1 = cell(1,num_timepoints);
         activity_matrices2 = cell(1,num_timepoints);

       for i=1:num_timepoints
           
            activity_matrices1{i} = obj.compute_pairwise_activity_timepoint(channel_names, channel_labels,i);
            activity_matrices2{i} = cseries2.compute_pairwise_activity_timepoint(channel_names,channel_labels, i);
%             
%            
%            minval = min(min(min(sensitivity_matrix1)),min(min(sensitivity_matrix2)));
%            maxval = max(max(max(sensitivity_matrix1)),max(max(sensitivity_matrix2)));
            
%            activity_matrices1{i} = obj.compute_pairwise_activity_timepoint(channel_names, channel_labels,i, 0, maxvals(i));
%            activity_matrices2{i} = cseries2.compute_pairwise_activity_timepoint(channel_names,channel_labels, i, 0, maxvals(i));
%            
       end
       
       
    

    end
    
    
    
    
    function [sensitivity_matrix] = compute_pairwise_activity_timepoint(obj, channel_names, channel_labels, time_point, varargin)
        
        figfile = sprintf('%s_time_%s_v2',obj.series_name,obj.datalabels{time_point});
        
        num_channels = length(channel_names);
        sensitivity_matrix = zeros(num_channels, num_channels);
       
        for i=1:num_channels
            for j=1:num_channels
                
                if(i==j) 
                    continue;
                end
                
                channel1_name = channel_names{i};
                channel2_name = channel_names{j};
                
                sensitivity_matrix(i,j) = obj.compute_activity(channel1_name, channel2_name, time_point,'mi_threshold',.01);
                
                    
            end
        end
        
%         figure
%         h = gca
%         
%         minval = min(min(sensitivity_matrix));
%         maxval = max(max(sensitivity_matrix));
%         
%         if(length(varargin)>0)
%             minval = varargin{1};
%         elseif(length(varargin)>1)
%             maxval = varargin{2};
%         end
%         
%         
%         colormap(pink);
%         
%         CLIM = [minval maxval];
%         imagesc(sensitivity_matrix, CLIM);
%         set(gca,'ytick',1:length(channel_labels));
%         set(gca,'yticklabel',channel_labels);
%         %xticklabel_rotate([1:length(channel_labels)],45,channel_labels);
%         saveas ( h, figfile, 'pdf' );
%         saveas ( h, figfile, 'fig' );
        
        
        
    end
       
    function [sensitivity_diff_matrix] = compute_pairwise_activity_difference_timepoint(obj, cseries2, channel_names, channel_labels, time_point, varargin)
        
        
        figfile = sprintf('%s_%s_compare_time_%s',obj.series_name, cseries2.series_name,obj.datalabels{time_point});
        num_channels = length(channel_names);
        sensitivity_diff_matrix = zeros(num_channels, num_channels);
       
        for i=1:num_channels
            for j=1:num_channels
                
                if(i==j) 
                    continue;
                end
                
                channel1_name = channel_names{i};
                channel2_name = channel_names{j};
                
                sensitivity_diff_matrix(i,j) = obj.compute_activity_difference(cseries2, channel1_name, channel2_name, time_point, 'mi_weight',.01);
                
                    
            end
        end
%         
%         figure
%         h = gca
%         
%         
%         
%         minval = min(min(sensitivity_diff_matrix));
%         maxval = max(max(sensitivity_diff_matrix));
%         
%         if(length(varargin)>0)
%             minval = varargin{1};
%         elseif(length(varargin)>1)
%             maxval = varargin{2};
%         end
%         
%         if(abs(minval)>abs(maxval))
%             maxval = -1*minval;
%         else
%             minval = -1*maxval;
%         end
%         
%         colormap(redgreencmap);
%         
%         CLIM = [minval maxval];
%         imagesc(sensitivity_diff_matrix, CLIM);
%         set(gca,'ytick',1:length(channel_names));
%         set(gca,'yticklabel',channel_labels);
%         %set(gca,'xaxisLocation','top')
%         xticklabel_rotate([1:length(channel_names)],90,channel_labels);
%        
%         saveas ( h, figfile, 'pdf' );
        
        
    end
    
    function [areas1, areas2] = pairwise_visualization_compare(obj, cseries2, channel1_name, channel2_name)
        
        
        
        
           num_datasets = length(obj.datasets);
        
           num_plot_rows = 2;
           num_plot_cols = num_datasets;
           maxy = 0;
           areas1 = zeros(1,num_datasets);
           areas2 = zeros(1,num_datasets);
         
           
          for j = 1:num_datasets
          

            
            cdataj = obj.datasets{j};
            cdatak = cseries2.datasets{j};
           
            channel1 = cdataj.name_channel_map(channel1_name);
            channel2 = cdataj.name_channel_map(channel2_name);
          
            subplot(num_plot_rows, num_plot_cols, j);    
            [points_x1, points_y1] = cdataj.pairwise_visualize(channel1_name,channel2_name, 'show_density');
            subplot(num_plot_rows, num_plot_cols, j+num_plot_cols);
            [points_x2, points_y2] = cdatak.pairwise_visualize(channel1_name,channel2_name, 'show_density');

            
            
            range_min = max(min(points_x1),min(points_x2));
            range_max = min(max(points_x1),max(points_x2));
   
            estimation_pts = linspace(range_min, range_max, 100);
            [points_y1_estimate] = interp1(points_x1, points_y1, estimation_pts, 'linear');
            [points_y2_estimate] = interp1(points_x2, points_y2, estimation_pts, 'linear');
   
            sensitivity = norm(points_y1_estimate-points_y2_estimate,1);
            areas1(j) = norm(points_y1_estimate,1);
            areas2(j) = norm(points_y2_estimate,1);
           
            
        
          end
    
         for j = 1:num_datasets
             
            subplot(num_plot_rows, num_plot_cols, j);
            xlabelstring = sprintf('%f',areas1(j));
            xlabel(xlabelstring);
             
            subplot(num_plot_rows, num_plot_cols, j+num_plot_cols);
            xlabelstring = sprintf('%f',areas2(j));
            xlabel(xlabelstring);
             
         end
       

    end
    
    
    function [node_difference, pt1, pt2] = compute_node_difference(obj, cseries2, channel_name, timepoint, varargin)
        
        options = 'percentage';
        if(length(varargin)>0)
            options = varargin{1};
        end
        
        %fold change 
        
        cdataj = obj.datasets{timepoint};
        cdatak = cseries2.datasets{timepoint};
        pt1 = cdataj.get_channel_mean(channel_name);
        pt2 = cdatak.get_channel_mean(channel_name);
        
        if(strcmp(options,'fold_change'))
            cdataj0 = obj.datasets{1};
            cdatak0 = cseries2.datasets{1};
        
         
             pt10 = cdataj0.get_channel_mean(channel_name);
    
             pt20 = cdatak0.get_channel_mean(channel_name);
             pt1 = pt1/pt10;
             pt2 = pt2/pt20;
        end
        
        node_difference = (pt1-pt2)/pt1;
        
        
    end
    
    function [edge_compare] = compute_edge_mis_compare(obj, cseries2, edges, timepoint)
    
      edge_entropys1 = zeros(length(edges),1);
      edge_entropys2 = zeros(length(edges),1);
      edge_compare = zeros(length(edges),1);
    
        for i=1:length(edges)
            
           [edge_entropys1(i),edge_entropys2(i)] = compute_centropy_compare(obj, cseries2, edge2{i,1}, edges{i,2}, timepoint);
          
           
        end
        
        edge_compare = edge_entropys1-edge_entropys2;
        
    end
    
    function [edge_pvals] = compute_edge_mis(obj, edges, timepoint)
    
      edge_pvals = zeros(length(edges),1);
      %have to fix this to compare fairly as well. 
      
    
        for i=1:length(edges)
            
           
            
            [edge_pvals(i)] = obj.compute_mi(edges{i,1},edges{i,2},timepoint);
           
        end
        
    end
    
    
   
   function [edge_entropys] = compute_edge_centropys(obj, edges, timepoint)
    
      edge_entropys = zeros(length(edges),1);
     
    
        for i=1:length(edges)
            
           
            
            [edge_entropys] = obj.compute_entropy(edges{i,1},edges{i,2},timepoint);
           
        end
        
    end
    
    
    
    function [edge_differentials, edge1, edge2] = compute_edge_differentials(obj, cseries2, edges, timepoint)
    
      edge_differentials = zeros(length(edges),1);
      edge1 = zeros(length(edges),1);
      edge2 = zeros(length(edges),1)
    
        for i=1:length(edges)
            
           
            
            [edge_differentials(i), edge1(i), edge2(i)] = compute_activity_difference(obj, cseries2, edges{i,1}, edges{i,2}, timepoint,'mi_threshold');
           %[edge_differentials(i), edge1(i), edge2(i)] = compute_activity_difference(obj, cseries2, edges{i,1}, edges{i,2}, timepoint);
           
            
        end
        
    end
    
    function [edge_differentials, edges1, edges2] = compute_all_edge_differentials(obj, cseries2, edges)
        
       edge_differentials = [];
       edges1=[];
       edges2=[];
       
       for i=1:length(obj.datasets)
           
          [new_differentials, newE1, newE2] = obj.compute_edge_differentials(cseries2, edges, i);  
          edge_differentials= [edge_differentials new_differentials];
          edges1 = [edges1 newE1];
          edges2 = [edges2 newE2];
       end
        
        
    end
    
    function [node_differentials, nodes1, nodes2] = compute_all_node_differentials(obj, cseries2, nodes)
        
       node_differentials = [];
       nodes1=[];
       nodes2=[];
       
       for i=1:length(obj.datasets)
           
          [new_differentials, newN1, newN2] = obj.compute_node_differentials(cseries2, nodes, i);  
          node_differentials= [node_differentials new_differentials];
          nodes1 = [nodes1 newN1];
          nodes2 = [nodes2 newN2];
       end
        
        
    end
    
    
    function [node_differentials, node1, node2] = compute_node_differentials(obj, cseries2, nodes, timepoint)
    
        node_differentials = zeros(length(nodes),1);
        node1 = zeros(length(nodes),1);
        node2 = zeros(length(nodes),1);
        for i=1:length(nodes)
            
            
           [node_differentials(i), node1(i), node2(i)] = obj.compute_node_difference(cseries2, nodes{i}, timepoint); 
           
        end
        
    end
    
   
    
      
   

    end
    

end
   
