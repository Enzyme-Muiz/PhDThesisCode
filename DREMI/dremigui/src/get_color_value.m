function [color_value] = get_color_value(value, range, colMAP)

        value = value/range;
        n = length(colMAP);
        scaled_value = ceil(value * n);
        if(scaled_value <= 1 ) 
            scaled_value = 1;
        end
        color_value = colMAP(scaled_value,:);
        
end

