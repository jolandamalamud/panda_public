function task = combine_sessions(D)
    
    ntrials = 96;
    session = [0 1 1];
    session_vector = repmat(session,96,1); 
    session_vector = session_vector(:);
    ids = str2num(vertcat(D.sjid));
    unique_ids = unique(ids);
    sess = str2num(vertcat(D.sess));
    task = struct();
    
    for s = 1:length(unique_ids)
        sub_sess = sess(ids == unique_ids(s));
        task(s).sess = session_vector;
        subject_data = D(ids == unique_ids(s));
        task(s).a = nan(ntrials*3, 1);
        task(s).s = nan(ntrials*3, 1);
        task(s).r = nan(ntrials*3, 1);
        for ss = 1:length(sub_sess)
            task(s).a((ss-1)*ntrials+1:ss*ntrials) = subject_data(ss).a;
            task(s).s((ss-1)*ntrials+1:ss*ntrials) = subject_data(ss).s;
            task(s).r((ss-1)*ntrials+1:ss*ntrials) = subject_data(ss).r;
        end
        task(s).sjid = num2str(unique_ids(s));
        task(s).Nch = sum(~isnan(task(s).a));
        task(s).sessions = sub_sess;
    end
    
%     for s = 1:length(unique_ids)
%         sub_sess = sess(ids == unique_ids(s));
%         task(s).data = D(ids == unique_ids(s));
%         for ss = 1:length(sub_sess)
%             task(s).data(ss).Nch = length(~(task(s).data(ss).a~=3));
%             task(s).data(ss).sess = repmat(session(ss), 96, 1);
%             task(s).data(ss).a(task(s).data(ss).a==3) = nan;
%             task(s).data(ss).r(task(s).data(ss).a==3) = nan;
%             task(s).data(ss).s(task(s).data(ss).a==3) = nan;
% %             task(s).data(ss).r(task(s).data(ss).r == worse(ss)) = -1;
% %             task(s).data(ss).r(task(s).data(ss).r == better(ss)) = 1;
%         end
%         task(s).sjid = num2str(unique_ids(s));
%         task(s).Nch = sum([task(s).data.Nch]);
%         task(s).sess = sub_sess;
%     end
        
end