classdef Bag
    properties
        elements
        last_idx
    end
    
    methods         
        function obj = Bag(cap)
           obj.elements = zeros(cap,1);
           obj.last_idx = 0;
        end
        
        function obj = add(obj,elem)
            obj.last_idx = obj.last_idx+1;
            obj.elements(obj.last_idx) = elem;
        end
    end
    
end

