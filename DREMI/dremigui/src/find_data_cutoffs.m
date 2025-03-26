function [minx, miny, maxx, maxy] = find_data_cutoffs(cdata, channel1_name, channel2_name, cutoff, num_slices)    

        channel1 = cdata.name_channel_map(channel1_name);
        channel2 = cdata.name_channel_map(channel2_name);
        
        X = cdata.data(:,channel1);
        Y = cdata.data(:,channel2);
        total_cells = size(X,1);
        
        
        [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices+1,[min(X) min(Y)],[max(X) max(Y)]);
   
        
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
   
   end