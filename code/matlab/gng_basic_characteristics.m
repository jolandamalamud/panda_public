function bc = gng_basic_characteristics(data)
    
    cr = [1 1 2 2]; % right response
    worse_outcome = [0 -1 0 -1]; % negative outcome
    
    for sj = 1:length(data) 
     
        to_late = data(sj).a == 3; ptolate(sj) = nanmean(to_late); 
        a = data(sj).a; a(to_late) = nan; 
        s = data(sj).s; 
        r = data(sj).r; r(to_late) = nan; 
        worse_outcome = sum(s == [2,4],2) * -1;
%         rt = data(sj).rt; rt(to_late) = nan;
        
        for ss=1:4
            i = s==ss;
            bc.rs(:,ss,sj) = r(i);
            bc.as(:, ss, sj)= a(i);
            bc.N(ss, sj) = length(a(i & ~isnan(a)));
            bc.pgo(ss, sj) = nanmean(a(i & ~isnan(a))==1);
            bc.pcorr(ss,sj) = nanmean(a(i & ~isnan(a))==cr(ss));
%             bc.rt(:, ss,sj) = rt(i);
%             bc.switch_incondition(:,ss,sj) = bc.as(1:end-1,ss,sj) ~= bc.as(2:end,ss,sj);
%             rs = r(i);
%             bc.switch_incondition_after_worse_last_trial(:,ss,sj) = bc.as(1:end-1,ss,sj) ~= bc.as(2:end,ss,sj) & r(find(i) - 1) == worse_outcome(find(i) - 1); %  & rs(1:end-1) == worse_outcome(ss) 
%             bc.switch_incondition_after_worse(:,ss,sj) = bc.as(1:end-1,ss,sj) ~= bc.as(2:end,ss,sj)& bc.as(1:end-1,ss,sj) ~= cr(ss); %  & rs(1:end-1) == worse_outcome(ss) 
%             bc.switch_incondition(:,ss,sj) = bc.as(1:end-1,ss,sj) ~= bc.as(2:end,ss,sj);
%             bc.stay_incondition(:,ss,sj) = bc.as(1:end-1,ss,sj) == bc.as(2:end,ss,sj);
%             rs = r(i);
%             bc.stay_incondition_after_better(:,ss,sj) = bc.as(1:end-1,ss,sj) == bc.as(2:end,ss,sj) & bc.as(1:end-1,ss,sj) == cr(ss); %& rs(1:end-1) == worse_outcome(ss)+1
%             
        end
                
%         beta = regress(rt, [ones(length(s),1), ~sum(s==[1,2], 2)]);
%         bc.mrt(sj) = beta(1);
        
        bc.pcorr_pav(1, sj) = bc.pcorr(1,sj) - bc.pcorr(3,sj);
        bc.pcorr_pav(2, sj) = bc.pcorr(4,sj) - bc.pcorr(2,sj);
        
    end
end