function p = ind_out_prob(s, a, r)
    
    % probabilities set up in task
    for k = 1:4; p(k,1) = 0.8; p(k,2) = 0.8; end
    
    
%     % individual probabilities of individual task runs
%     rew = [1 0 1 0]; pun = [0 -1 0 -1];
%     cr = [1 1 2 2]; wr = [2 2 1 1];
%     for k = 1:4
%         if sum(a(s==k)'== cr(k)) > 3
%             p(k,1) = sum(r(s==k) == rew(k) & a(s==k) == cr(k)) / sum(a(s==k) == cr(k));
%         else
%             p(k,1) = 0.8;
%         end
%         if sum(a(s==k)'== wr(k)) > 3
%             p(k,2) = sum(r(s==k) == pun(k) & a(s==k) == wr(k)) / sum(a(s==k) == wr(k));
%         else
%             p(k,2) = 0.8;
%         end
%     end
end