addpath(genpath('/analyse/Project0234/angus/colldiag/'))
clear all;clc;

% A slow event-related design - design matrices will have variable total volumes

% subject number
sub_num = 2;

time = 300; % max scanning time in seconds
tr = 2; % sample rate in seconds
t = 1:tr:time; % measurements
h = gampdf(t,6) + -.5*gampdf(t,10); % hrf model
h = h/max(h); % scale hrf to max amplitude 1

start_base = 6/tr; %samples, baseline before first onset
end_base = 12/tr; %samples, baseline after last onset
conditions = 7;
reps = 8;
runs = 4;
colThr = 0.3;

r = unifrnd(5,9,[1,1000]); % for choosing ISI
r_count = 1;

while r_count <= runs % get 4 designs below thresh

        % Produce the pseudo-randomised stimulation presentation order

        % Generate a vector of trial order.
        trials = repmat(1:conditions,reps,1); 
        trials = trials(:); % vectorise
        trials = trials(randperm(length(trials))); % randomise the order
        order = [];
        tmp = datasample(trials,reps*conditions,'Replace',false); %randomly sample without replacement
        order = tmp; % the new order

        % Generate some design matrix.
        A = 0;

        % start / start again
        count = start_base;
        design = zeros(time/tr,1);

        %first trial
        count = count+1;
        design(count,1) = 1;
        stim_c = 1;

        while stim_c < reps*conditions

            % choose an ISI
            isi = datasample(r,1);
            isi_out(stim_c,r_count) = isi;
            count = count + floor(isi/tr) + 1;
            design(count,1) = 1;
            stim_c = stim_c + 1;

        end

        out = design;
        A = 0; design_matrix = [];
        
%         if size(design,1) == time/tr

            % create design matrix and convolve with hrf model
            [A] = find(order == 1); 
            [B] = find(order == 2); 
            [C] = find(order == 3); 
            [D] = find(order == 4);
            [F] = find(order == 5);
            [G] = find(order == 6);
            [H] = find(order == 7);
            [E] = find(out == 1);
            DM = zeros(time/tr,2);
            DM(E(A),1) = 1;
            DM(E(B),2) = 1;
            DM(E(C),3) = 1;
            DM(E(D),4) = 1;
            DM(E(F),5) = 1;
            DM(E(G),6) = 1;
            DM(E(H),7) = 1;
            design_matrix(:,1) = conv(DM(:,1),h);
            design_matrix(:,2) = conv(DM(:,2),h);
            design_matrix(:,3) = conv(DM(:,3),h);
            design_matrix(:,4) = conv(DM(:,4),h);
            design_matrix(:,5) = conv(DM(:,5),h);
            design_matrix(:,6) = conv(DM(:,6),h);
            design_matrix(:,7) = conv(DM(:,7),h);
            design_matrix(size(DM,1)+end_base:end,:) = [];
            design_matrix(:,end)=1; % replace null condition with run constant
            
            [T] = evalc('[sValue,condIdx,VarDecomp] = collintest(design_matrix);');
            Vd = max(max(VarDecomp(1:6,1:6)));
            disp(Vd)

            if Vd <= colThr
                store(r_count).dm = design_matrix;
                store(r_count).seq = order;
                store(r_count).design = design;
                store(r_count).isi = isi_out(:,r_count);
                disp(sprintf('Found design %d',r_count))
                r_count = r_count+1;
            end
%         else
%         end
end

%% Output files

if ~exist(sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/',sub_num))
    mkdir(sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/',sub_num));
    fileattrib(sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/',sub_num),'+w','a');
end

for rr = 1:runs
    outname=sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/RUN%d_SEQUENCE.txt',sub_num,rr);
    fid=fopen(outname,'w');
    
    % sequence
    sequence = store(rr).seq;
    
    for jj = 1:size(sequence,1)
         fprintf(fid, sprintf('%d',sequence(jj,1)));
         if jj ~= size(sequence,1)
             fprintf(fid, ',');
         end
    end
    fprintf(fid, '\n');
    fclose(fid);

    % jitter
    outname=sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/RUN%d_JITTER.txt',sub_num,rr);
    fid=fopen(outname,'w');
    
    jitter = store(rr).isi;
    for jj = 1:size(jitter,1)
         fprintf(fid, sprintf('%d',jitter(jj,1)));
         if jj ~= size(jitter,1)
             fprintf(fid, ',');
         end
    end
    fprintf(fid, '\n');
    fclose(fid);
    
    % design
    outname=sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/RUN%d_DESIGN.txt',sub_num,rr);
    fid=fopen(outname,'w');
    
    design = store(rr).design;
    for jj = 1:size(design,1)
         fprintf(fid, sprintf('%d',design(jj,1)));
         if jj ~= size(design,1)
             fprintf(fid, ',');
         end
    end
    fprintf(fid, '\n');
    fclose(fid);
    
end

save(sprintf('/analyse/Project0385/angus/DESIGNS/SUB%d/subject%d.mat',sub_num,sub_num),'store')

%% 
% Vd = 1;
% while Vd > 0.3
%     sb_matrix = zeros(size(design_matrix,1),reps*conditions);
%     mm = start_base;
%     for ii = 1:reps*conditions
%         mm = mm + floor(isi_out(ii)/tr);
%         tmp = zeros(size(design_matrix,1),1);
%         tmp(mm) = 1;
%         tmp = conv(tmp,h);
%         sb_matrix(:,ii) = tmp(1:size(design_matrix,1));
%         mm = mm+2;
%     end
%     sb_matrix(:,end+1) = 1;
% 
% 
%     [T] = evalc('[sValue,condIdx,VarDecomp] = collintest(sb_matrix);');
%     Vd = max(max(VarDecomp(1:reps*conditions,1:reps*conditions)));
%     disp(Vd)
% end