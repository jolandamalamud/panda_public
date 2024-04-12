function [a,r] = generatera(pa, s, p);

    c = 1;
    ago = 1; 

    a = sum(rand>cumsum([0 pa]));

    if s==1	% go to win 
        if a == ago; 	r = c*(rand<p(1,1));
        else; 			r = c*(rand>p(1,2)); 
        end
    elseif s==2; 	% go to avoid 
        if a == ago; 	r = -c*(rand>p(2,1));
        else; 			r = -c*(rand<p(2,2));
        end
    elseif s==3		% nogo to win 
        if a == ago; 	r = c*(rand>p(3,2)); 
        else; 			r = c*(rand<p(3,1));
        end
    elseif s==4; 	% nogo to avoid 
        if a == ago; 	r = -c*(rand<p(4,2));
        else; 			r = -c*(rand>p(4,1));
        end
    end
    
%     if s==1	% go to win 
%         if a == ago; 	r = c*(rand<p(1,1));
%         else 			r = 0; 
%         end
%     elseif s==2; 	% go to avoid 
%         if a == ago; 	r = 0;
%         else 			r = -c*(rand<p(2,2));
%         end
%     elseif s==3		% nogo to win 
%         if a == ago; 	r = 0; 
%         else 			r = c*(rand<p(3,1));
%         end
%     elseif s==4; 	% nogo to avoid 
%         if a == ago; 	r = -c*(rand<p(4,2));
%         else 			r = 0;
%         end
%     end
    
end