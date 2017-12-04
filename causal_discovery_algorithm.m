% Causal discovery algorithm for quantum systems
% Requires: Functions from Tony Cubitt: syspermute, tensor, TrX
%           My functions: traceout, traceout2, traceoutcom, traceoutcomN
%           Other functions: cprintf
% Author: Christina Giarmatzi
% See Manual for instructions.


clear all

% Input of the user.
tag = '40'; % the name tag of dmat[tag], W[tag] and subdim[tag].
folder = '/Users/...';% The folder that has the input

% Imports
Wname = ['W' tag '.dat'];
dmat = ['dmat' tag '.dat'];
dmator = load([folder dmat]);% To record the original matrix dmat. Later dmat might change.
dmat = dmator;
N = size(dmat,1);% The number of parties.
W = (load([folder Wname]));
sname = ['subdim' tag '.m'];
subdim{N} = [];
run(sname);
subdimor = subdim;% To record the original subdim. It might change later.
epsilon = 10e-5;

% Close any open biograph viewer from previous test.
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))

tic % Start the clock.

spart = []; % The parties that have subsystems.
for xx = find(~cellfun(@isempty, subdim));
spart = [spart xx];
end

% Indexing will keep track of the right labels of the subsystems. If subs 1 is open and traced out, then subs 2 will be involved in a link and 
% indexing makes sure that it will be called subs 2 although in the rest of the code it will be thought as subs 1.
indexing{N} = []; % in case there are no subsystems and indexing never gets defined
for uu = spart;
indexing{uu}= zeros(1,length(subdim{uu}));% for the parties that have subsystems, indexing is defined
end

dim = reshape(dmat',1,2*N); % A row of [di1 do1 di2 do2 ...]
dimor = dim;% dim original, because dim will change

% One dimensional systems: If it is an input: insert identity/2. If it's an output: insert identity. And replace the 1s with 2s in dim.
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));
for onedim = find(dim==1);
    if mod(onedim,2);% if the one dimensional system is an input
        party = (onedim+1)/2;
        W = syspermute(tensor(eye(2)/2, W), [2 1 3], [2 prod(dim(1:2*(party-1))) prod(dim(2*party:length(dim)))]);
        dim(onedim) = 2;
    else% if the one dimensional system is an output
        party = onedim/2;
        W = syspermute(tensor(eye(2), W), [2 1 3], [2 prod(dim(1:2*party-1)) prod(dim(2*(party+1)-1:length(dim)))]);
        dim(onedim) = 2;
    end
end
dmat(dmat==1) = 2; % Fix also dmat: replace the one dimensions with two.

Norm = prod(dmat(:,2));% Product of dimensions of all input systems.
Trace = trace(W);


% Normalization of the process matrix. 
if abs(Norm-real(Trace))>10e-4;
     cprintf('red', 'Process matrix got normalized');
     W = (Norm/real(Trace))*W;
end

% Defining useful matrices.
d_all = prod(dim); % dimension of W(d_all,d_all)
NLast = zeros(d_all, d_all, N); % Nlast(:,:,N) is a process matrix compatible with Nth party to be last.
Sets = zeros(N,N); % In the first row there will be the (numbering of) parties that are 1st Last, in the second will be the 2nd Last end so on.
W_rest = zeros(d_all, d_all, N); % the W after tracing out the parties in Nth set
parties_sum = 0;% Counts the parties used up so far, after each set
Wtest2 = zeros(d_all, d_all, N); % No use really. (The N numbering is tricky.) The process matrices after tracing out the nth and all the previous parties. 
% In the ordering that the code process them.
d_set = zeros(1,N); % Will have the dimensions of in-out systems of the parties in set=1 set=2 etc [dset1 dset2 ... ]

term{N,N+1} = []; %For every arrow/last party a term will be created to check Markovianity of W.
subs = sum(cellfun(@sum, subdim));% The sum of all the dimensions of subsystems.
max_sub = max(cellfun(@length, subdim));% The number of max subsystems composing one output system
tot_subs = sum(cellfun(@length, subdim)); % The number of all subsystems.
n_subs = length(find(~cellfun(@isempty, subdim))); % The number of output systems that are divided to subsystems.
counter1 = zeros(1,N);% The position in Sets from which parties start occupy. 
remaining_parties = [1:N];
if max_sub;% If the total number of systems is not zero, subterm{} will be created.
    subterm{N,N,max_sub} = []; % For every arrow for a subsystem, a term will be created.
end


%***************************************************************************************************************************************************
%******************************************* Causal order ****************************************************************************************
%***************************************************************************************************************************************************


Wtest1 = W; % just a buffer 
for set = 1:N; % the consecutive sets. set = 1 is the group of parties that are last.
    counter = counter1(set); % Pointer for the next position available in Sets(set, counter) to write down the number of the party who is in the consecutive 'set'.

    for n = remaining_parties; % Check all parties if they belong to set (without all the used ones)
        NLast(:,:,n) = traceout(Wtest1, 2*n, dim);% The traceout function normalizes the output matrix.
        if abs(Wtest1 - NLast(:,:,n))< 10e-4; % Asks the question: did the Wtest1 already had identity on that output / Is party n last?
        counter = counter +1;% Increase the position of the pointer in the row 'set'.
        Sets(set, counter) = n; % This matrix has in the first row all the parties that are last (set 1), in the second is the set 2 and so on.
        remaining_parties(remaining_parties==n) = [];% Remove n party from the remaining parties.
        end
    end
    % Basically if Sets(set,:)==zeros(size(Sets,1)) then it's non causally ordered.
    
    % Here the parties in 'set' have been established. 
    % Counter is the number of parties in the current set.
    d_set(set) = prod(dim(2*Sets(set,1:counter)-1))*prod(dim(2*Sets(set,1:counter))); % The product of dimensions of in-out systems of the parties in the 'set'.
    the_parties{set} = Sets(set,1:counter); % A cell array with the sets and the parties. length(the_parties{set}) is the number of parties in the set. 
    parties_sum = parties_sum + counter; % # parties used up so far.
    for m = the_parties{set}; % m = only the parties in the 'set'.
        Wtest2(:,:,m) = traceout2(Wtest1, [2*m-1 2*m], dim)*dim(2*m-1); %traceout normalizes by the dimension of the system that is traced out.
        % When tracing out a party, we need to renormalize the W by dividing with the dimension of the output system. That's why I'm *dim(input system).
        Wtest1 = Wtest2(:,:,m); % the new W is the one after tracing out n.
    end
    if parties_sum == N; % If # parties processed reaches N.
        sets = set; % The total number of sets is the 'set' in which the loop stopped.
        break % Stop the loop.
    end
end
if parties_sum~=N;
    cprintf('red', 'Non-causally ordered process matrix');
    return
end
if sets==1;
    %Then all parties are independend the parties in set=1 will have been written twice
    the_parties{1} = unique(the_parties{1});
end
% Remove zeros in the sets to present it on command window
the_sets = Sets;
the_sets(sets+1:size(Sets,1),:) = []; % delete the lines with zeros: sets is the number of sets, below that line Sets will have zeros
max_set = max(cellfun(@length, the_parties)); % this is the max parties that belong to a set. Sets's rows will be occupied up to that row.
the_sets(:,max_set+1:size(Sets,2)) = [];
the_sets
disp(['Time ' num2str(toc)])
fprintf('\n');


%**********************************************************************************************************************************
%*************************************************** Primal causal arrows  **********************************************
%**********************************************************************************************************************************

% Primal arrows are the arrows between parties that belong to adjacent sets (1-2, 2-3, etc)
% First look for open systems (identity on the output subsystem of a party) and trace them out of the W 
% Then causal arrows are established, between either output systems or output subsystems of party n2 and input system of party n1.
% Trace out input of n1, look for identity for the output of n2.
% Formally this is : Tr_Ai - Tr_AiBo < tolerance 
% Note: the code does not distinguish between two arrows "A->D and C>D" _and_ "(A and C)->D". In both cases it says A->D and C->D.
tic
dim_remain = dim;% Whenever a subsystem is being used in a causal arrow, dim_remain(i) = dim_remain(i)/subdim{i}(j). When a system is used dim_remain = 1.
% In the end, dim_remain shoud be a matrix ones(1,2*N)
% primal_arrows = zeros(N*(N-1)/2, 2);% Each row contains two parties whith a causal arrow from first to the second (the max amount of links between N dots is N*(N-1)/2).
tot = length(find(cellfun(@isempty, subdim)));
%tot is the number of parties with no output subsystems 
tot_rem = sum(cellfun(@(x) numel(x),subdim)) + tot; % The number of subsystems plus the number of output systems with no subsystem
poss_arr = tot_rem*(tot_rem-1)/2; % The number of possible arrows, given that all the possible ways of connecting n dots is n*(n-1)/2.
primal_arrows = zeros(poss_arr*(poss_arr-1)/2, 2);
store = zeros(tot_rem,3);% Here the arrows will be stored that correspond to output subsystems [n2 n1 p] the output party, input party and the index of the subsystem used.
storeoriginal = store;
subdim_original = subdim; % subdim_original must remain intact for the right labeling of subsystems
% Sanity check for the number of parties in each set (and the things used).
% 2 ways to calculate the # of parties in each set.
parties = zeros(sets,1);
for set = 1:sets;
parties(set) = length(the_parties{set});
end
Sets_logic = Sets~= 0; % A logic matrix with ones indicating a party on the set and zeros that have remained from the definition of Lasts
parties2 = sum(Sets_logic,2);
parties2(parties2==0) = [];
parties(parties==0) = [];
if parties2~=parties;
    disp('error in the calculation of the number of parties for each set');
end
insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end));
k = 0;% counter of arrows
counter2 = 0;
subopen{N}= [];
saveindex{N} = [];

for h = 1:N;% Only for output systems (sys2).
    subopen{h} = zeros(1, length(subdim{h})); % for every non-empty subdim{h}, create a zero cell array of the same length
end
subdim_remain = subdim;% Its elements will be removed as soon as they get used.
subdim_remain0 = subdim_remain; % Its elements will be zero as soon as they get used.
Setin = Sets;
Setout = Sets;
for set1 = 1:sets-1;% Checks sets: 1-2, 2-3, 3-4. for 4 sets
    for set2 = set1 + 1;
        for n1 = Setin(set1, 1:parties(set1)); % The parties inside set1, looking at their inputs at position dim(2*n1-1).
            sys1 = 2*n1-1; 
            for n2 = Setout(set2, 1:parties(set2)); % The parties inside set2, looking at their outputs  at position dim(2*n2).
                if n2>0;
                sys2 = 2*n2; 
                
                % The following constructs dim_prov = dim with inserted subsystem's dimensions in the place of the output of n2.
                % To be used in the next for_loops to check for causal arrows.
                
                dim_prov = dim;
                    if ~isempty(subdim_remain{n2});% If n2 has remaining subsystems
                        
                        ls2 = length(subdim{n2})-1;
                        dim_prov(sys2) = [];% Remove the sys2 element of dim_prov
                        dim_prov_new{1} = insert(subdim{n2}, dim_prov, sys2-1);%This is the right list of dimensions when checking for causal arrows
                        
% ******************************************************** Open output subsystems and remove them from W *******************************************************

                        for pout = 0:ls2;% Pointer to check all the subsystems.
                            openin = traceout(W, sys2+pout, dim_prov_new{1}) - W;% Zeros if there is an open output subsystem.
                            if abs(openin - zeros(d_all, d_all)) < epsilon;% If there is an open subsystem.
                              subopen{n2}(pout+1) = subdim{n2}(pout+1);% Store its index and dimension from subdim.
                              fprintf('There are open subsystems: %d of party %d of dimension %d', pout+1, n2, subdim{n2}(pout+1));
                              cprintf('green', '  --- *** ---');
                              fprintf('\n');
                              indexing{n2}(pout+1) = subdim{n2}(pout+1);% Keep the original labelling and dimension of the open subsystems.                              
                              dim_remain(sys2) = dim_remain(sys2)/subdim{n2}(pout+1); % Update the remaining dimensions removing the dim of open subsystem
                            end

                        end
                        
                        % subopen{n2} has non-zeros when there is an open subsystem with the number being the dimension
                        o = find(subopen{n2},1);% subopen has dimensions of the open subsystem. o is the first element: the first open subsystem
                        % The following will remove the subsystem from W, updating all the info about the dimensions.
                        % Finally the system is replaced by a 0 in the subopen, and the procedure starts again o = find(subopen{n2},1).
                        count0 = 0;
                     while ~isempty(o);
                        % the dimension of this open subsystem is removed. To keep track of the right index, we use count0.
%                         The dimensions are removed from left to right, so each time the new index will be the current -1.
                              o = o-count0;
                              subdim{n2}(o) = [];%  Remove it from subdim{n2}.
                              subdim_remain{n2}(o) = [];%  Remove it from subdim_remain{n2} too.
                              subdim_remain0{n2}(o+count0) = 0;
                              W = TrX(W, sys2-1+o, dim_prov_new{1})/dim_prov_new{1}(sys2-1+o);% Remove it from W and renormalize
                              if dim_prov_new{1}(sys2-1+o)~= indexing{n2}(o+count0); % Sanity check.
                                  cprintf('red', 'Sanity check failed for adressing the open subsystems.');
                              end
                              d_all = d_all/dim_prov_new{1}(sys2-1+o);% Renormalize d_all
                              dim(sys2) = dim(sys2)/dim_prov_new{1}(sys2-1+o);% Renormalize dim
                              dim_prov_new{1} = insert(subdim{n2}, dim_prov, sys2-1);% Create the new dim_prov_new with the updated subdim
                              subopen{n2}(o+count0) = 0;%Remove the information about the open subsystem that has just been removed to avoid falling into the for loop again.
                              o = find(subopen{n2},1);
                              count0 = count0+1;

                     end
                    else % if there are _no_ subsystems in sys2
                        dim_prov_new{1} = zeros(1,length(dim_prov));   
                    end
                       if dim_prov_new{1} ~= zeros(1, length(dim_prov_new{1}));
                           dim_prov = dim_prov_new{1};% if it is zeros, dim_prov will remain the original.
                       end
                       if ismember(n2,spart);
                            if ~isempty(indexing{n2}) & indexing{n2}==zeros(1,length(indexing{n2}));%Empty the zero matrices
                                indexing{n2}=[];
                            end
                       end

% ********************************************************* Check for primal arrows **********************************************
                                
            if ~isempty(subdim_remain{n2});  % If there _are_ subsystems in n2 (after having removed the open and the previously found arrows) 
                ls22 = length(subdim{n2})-1; % this is for the amount of positions to be moved 
                % Before when it was ls2 :I am checking all the subsystems whether some have been used or not
                ww=0;
                
                for pout = find(subdim_remain{n2})-1;
                    ww=pout+1;
                    if ~isempty(indexing{n2})
                      o2 = find(indexing{n2}==0);% For the correct labelling of the subsystems, 
% Here, indexing{n2} has all zeros apart from the places where there was an open subs, where it has the dim of that subs.
                    else
                      o2(ww) = pout+1;%if there were no open systems
                    end
                      if ((n1<n2) & abs(traceout(W, sys1, dim_prov) - traceout2(W, [sys1 sys2+pout], dim_prov) < epsilon)) | ...
                      ((n2<n1) & abs(traceout(W, sys1+ls22, dim_prov) - traceout2(W, [sys1+ls22 sys2+pout], dim_prov) < epsilon));
                  savep = [pout sys1 sys2 ls2 dim_prov];
                  
                  fprintf('Link from subsystem %d of party %d to party %d.',o2(ww), n2, n1);
                  cprintf('cyan', '  --- *** ---');
                  fprintf('\n');
                  k = k+1;
                  counter2 = counter2 +1;
                  store(counter2,:) = [n2 n1 pout+1];% NOT YET CLEAR IF IT'S CORRECT
                  storeoriginal(counter2,:) = [n2 n1 o2(ww)];
                  primal_arrows(k,:) = [n2 n1];% Store the info into a matrix [n1 n2 ; n1 n2; ... ]. Arrow from n2 to n1.
                  dim_remain(sys2) = dim_remain(sys2)/subdim{n2}(pout+1);
                  subdim_remain0{n2}(o2(ww)) = 0;% This one has zero on the used and open systems.
                  subdim_remain{n2}(pout+1) = 0;% This one has the open systems removed and zero on the used systems.
               if subdim_remain{n2} == zeros(1,length(subdim_remain{n2}));
                   Setout(set2,Setout(set2,:)==n2)=0; 
               end
% saveindex{n2} is the subsystems of n2 that have been used up in an arrow.
                  saveindex{n2} = [saveindex{n2} o2(ww)];% it needs to have the original index to be used with indexing{}
%  *******************  Markovianity    *******************
% Here traceoutcom does not do the right job with normalization. This is because to check whether an
% ouput system is being traced (to add to the normalization) it checks if every system that is traced out
% is an even number. But this does not work well with dim_prov that has extra systems inserted,
% which mess up the numbering of input output systems from odd even.
% Therefore, traceoutcomN is used that does not normalize the output.
% Then normalize the term by /by its trace and multiplying with its dim(output).
% link is from sysbsystem of party n2 to party n1.
                  if n1<n2;
                      buff1 = traceoutcomN(W, [sys1 sys2+pout], dim_prov);% Unnormalized term.
                      buff2 = buff1/trace(buff1)*dim_prov(sys2+pout);% Normalized term.
                      subterm{n2,n1,o2(ww)} = syspermute(buff2, [2 1], [dim_prov(sys1) dim_prov(sys2+pout)]);% Swap it to make it output-input.
                  else 
                      buff1 = traceoutcomN(W, [sys1+ls22 sys2+pout], dim_prov);% Unnormalized term
                      subterm{n2,n1,o2(ww)} = buff1/trace(buff1)*dim_prov(sys2+pout);% Normalized.
                  end
                  
                      end
                end
% 
               
            else% if there are _no_ subsystems in n2
               
                if ( abs(traceout(W, sys1, dim_prov) - traceout2(W, [sys1 sys2], dim_prov)) < ( epsilon ) ) ;
                  fprintf('Link from party %d to party %d.' , n2, n1);
                  cprintf('blue', '  --- *** ---');
                  fprintf('\n');
                  k = k+1;
                  counter2 = counter2+1;
                  primal_arrows(k,:) = [n2 n1];% Store the info into a matrix [n1 n2 ; n1 n2; ... ]. Arrow from n2 to n1.
                  store(counter2,:) = [n2 n1 0];
                  storeoriginal(counter2,:) = [n2 n1 0];
                  dim_remain(sys2) = 1;
                  Setout(set2,Setout(set2,:)==n2)=0;
                  % Markovianity
                  term{n2,n1} = traceoutcom(W, [sys1 sys2], dim_prov);
                  if n2>n1;
                      term{n2,n1} = syspermute(term{n2,n1}, [2 1], [dim(sys1) dim(sys2)]);
                  end
                
                end
            end  
                end
            end
        end
    end
end

%%%% The following is not necessary. In fact it produces errors. The subdim_remain keeps track of the used systems and then the code only 
% checks for arrows on the output systems of parties that have a nonzero entry in the subdim_remain!
%for n2=1:N; % it's done again after the secondary arrows
 %   if ~isempty(indexing{n2});
 %    indexing{n2}(saveindex{n2}) = 1;% Now the indexing will have 0s only when the systems have not only been not opened but also not used.
 %   end
% Here, k is the number of primal arrows detected by the code.
% Need to remove the zeros from the primal_arrows

primal_arrows(k+1:length(primal_arrows),:) = [];
fprintf('\n')
disp([num2str(k) ' primal arrows']);
primal_arrows
disp(['Time ' num2str(toc)])
fprintf('\n');

%***********************************************************************************************************
%********************************* Secondary causal arrows *************************************************
%***********************************************************************************************************

%  Same technique as for the primal arrows is applied, except that now the code is looking for links between a party and 
% _all_ the parties that belong the sets 'lower' than the set in which A belongs to.

tic
subdim_remain1 = subdim_remain; % to update it every time there is an arrow, to remove the subsystem involved
% ouputs the links with the _right_ labels of subsystems. 
% Now, dim_remain has zeros at the place where the output systems were used up for causal arrows.
% Also, subdim_remain{n2} has zeros in the place where the output subsystems were used up.
% Remove all zeros from subdim_remain, so that it has only remaining subsystems, for the next code to work.
% (maybe there's a better way to write this, with a cellfun)

for i = 1:N;
    if ~isempty(subdim_remain{i}) & (subdim_remain{i} == zeros(1,length(subdim_remain{i})));
        % If it's not empty and has zeros, empty it.
        subdim_remain{i} = [];
    end
end
% Now indexing has the dimension of output systems on the open ones, and 1 on the used up ones.
% So if there are no zeros remaining on indexing, then no output arrow can be found
dim_remain_out = dim_remain(2:2:2*N); % An array with only the output remaining subsystems.
logic_remain_out = dim_remain_out ~= 1 ; % A logical array for the output ones.
party_out = find(logic_remain_out);% the parties that have remaining output systems
tot_rem = sum(cellfun(@(x) length(x),subdim_remain)); % The number of subsystems remaining. Also works with numel(x) (number of elements)
if tot_rem;
    num_sec_arr = 1;
else
    num_sec_arr = tot_rem*(tot_rem-1)/2;
end
secondary_arrows = zeros(num_sec_arr, 2);
subdim_remained_after_primal = subdim_remain; % Just for record if the right index of non-zero remaining elements.
%  Links from n2 to n1
l = 0; % counter
for set1 = 1:sets-2; % each set
    for set2 = set1+2:sets; % each remaining set but not the immediate next
        for n1 = the_parties{set1}; % The parties from set1
            sys1 = 2*n1-1;
            for n2 = the_parties{set2}(find(ismember(the_parties{set2}, party_out))); % party_out has the parties with unused output systems or subsystems.
                sys2 = 2*n2;
                if dim(sys2)~=1;
                parties = the_parties{set2}(find(ismember(the_parties{set2}, party_out)));
                
 % The following constructs dim_prov = dim with inserted subsystem's dimensions in the place of the output of n2.
 % The subdim will be inserted, with all the original dimensions (even of used up systems). This is because ...
 % when tracing out these subsystems from W, these used up subsystems were not traced out from W.
 
                dim_prov = dim;% If there are no subsystems, the following 'if' will not be executed.
                
% ******************************************** Check for secondary causal arrows ******************************************************

            if ~isempty(subdim_remain{n2});  % If there _are_ still subsystems in sys2  
                dim_prov(sys2) = [];% Remove the sys2 element of dim_tes
                dim_prov_new{2} = insert(subdim{n2}, dim_prov, sys2-1);
                dim_prov = dim_prov_new{2};
                ls2 = length(subdim{n2})-1;%the 0 is because it has 0 when the subsystem has been used
%                 RM subdim_remain{n2} has 0 when the system was used up and no open subsystems
              w = 0;
              
                for pout = find(subdim_remain{n2})-1; % look only for the remaining subs
                     w = pout+1;
                    if ~isempty(indexing{n2});
                      o4 = find(indexing{n2}==0);% Zero when the subsystem was not open
                    else
                      o4(w) = pout+1;
                    end
                      if ((sys1<sys2) & (abs(traceout(W, sys1, dim_prov) - traceout2(W, [sys1 sys2+pout], dim_prov)) < epsilon)) | ...
                      ((sys2<sys1) & (abs(traceout(W, sys1+ls2, dim_prov) - traceout2(W, [sys1+ls2 sys2+pout], dim_prov)) < epsilon));
                  
                 savep2 = pout;
                  fprintf('Link from subsystem %d of party %d to party %d.', o4(w), n2, n1);
                  cprintf('cyan', '  --- *** ---');
                  fprintf('\n');
                  l = l+1;% counting the number of arrows to display later
                  counter2 = counter2+1;
                  store(counter2,:) = [n2 n1 pout+1];
                  storeoriginal(counter2,:) = [n2 n1 o4(w)];
                  secondary_arrows(l,:) = [n2 n1];% Store the info into a matrix [n1 n2 ; n1 n2; ... ]. Arrow from n2 to n1.
                  dim_remain(sys2) = dim_remain(sys2)/subdim_remain{n2}(pout+1);% Update the dim_remain
                  subdim_remain0{n2}(o4(w)) = 0;% This one has zero on the used and _open_ systems.
                  subdim_remain{n2}(pout+1) = 0;% This one has the open systems removed and zero on the used systems.
                  saveindex{n2} = [saveindex{n2} o4(w)];
% NEW ADDITION-from here
                  dim_remain_out = dim_remain(2:2:2*N); % An array with only the output remaining subsystems.
                  logic_remain_out = dim_remain_out ~= 1 ; % A logical array for the output ones.
                  party_out = find(logic_remain_out);% the parties that have remaining output systems
% NEW ADDITION-to here  

% Markovianity
                  if n1<n2;
                      buff1 = traceoutcomN(W, [sys1 sys2+pout], dim_prov);% Unnormalized term.
                      buff2 = buff1/trace(buff1)*dim_prov(sys2+pout);% Normalized term.
                      subterm{n2,n1,o4(w)} = syspermute(buff2, [2 1], [dim_prov(sys1) dim_prov(sys2+pout)]);
                  else
                      buff1 = traceoutcomN(W, [sys1+ls2 sys2+pout], dim_prov);% Unnormalized.
                      subterm{n2,n1,o4(w)} = buff1/trace(buff1)*dim_prov(sys2+pout);% Normalized.
                  end
                      end
                end

            else% if there are _no_ subsystems in sys2, this is not possible. Distant arrows are always secondary.
                cprintf('red','Error: Open output system of %d was detected for a secondary arrow', n2);
                fprintf('\n');
                break
            end                   
                end
            end
        end
    end
end

for n2=spart;
 if ~isempty(indexing{n2});
     indexing{n2}(saveindex{n2}) = 1;% Now the indexing will have 0s only when the systems have not only been not open but also not used.
 end
end


fprintf('\n')
disp([num2str(l) ' secondary arrows']);% l is the number of secondary arrows.
if ~isempty(secondary_arrows);
    for iii = size(secondary_arrows,1);
        if secondary_arrows(iii,:)== [0 0];
            secondary_arrows(iii,:)=[];
        end
    end     
end
secondary_arrows
%

%******************************** Unused subsystems and causally independent parties *******************************************

for i = 1:N;
    if ~isempty(subdim_remain{i}) & (subdim_remain{i} == zeros(1,length(subdim_remain{i})));
        % If it's not empty and has zeros, empty it.
        subdim_remain{i} = [];
    elseif ~isempty(subdim_remain{i}) && ismember(i,the_parties{1}); %if there are unused systems from last parties, it makes sense 
        subdim_remain{i} = [];
        %for the Markovianity step
    end
end

% Remember, indexing{N}=[] and indexing{spart}=zeros(1,length(spart)) were initially defined.
% Then indexing{spart}(i) has the dimension of the open output systems which is never 1.
% Then it has ones, when the system was used in a causal arrow.
% Therfore, indexing{i} should never have any zeros, apart from the last parties whose outputs are never used.
for n2 = spart;% the parties that have subsystems.
        if ~ismember(n2, the_parties{1}) & ~isempty(find(indexing{n2}==0,1));% if n2 is not last and there is one or more zeros in indexing{n2}
            cprintf('red','There are unused output subsystems of party %d.', n2)
            fprintf('\n')
        end
end
% No need to check for unused input systems as it is fine if there are. 

% Now check for unused input systems and if they are last, they are causally independent.
% dim_remain = dim except that it has 1 in the output systems that were used.
% Now put 1 in the output systems of the last parties.
for op = the_parties{1};
    dim_remain(2*op) = 1;
end
% Now every input system that has been used by a causal arrow once, put its dimension as 1.
total_arrows = [primal_arrows;secondary_arrows];
ti = total_arrows(:,2);% the parties that have received an arrow. But there are double entries here.
[dummy, I] = unique(ti,'first');%Cool way to remove the sorting that unique does.
u = ti(sort(I));
usedinputs{1} = u';% These are the parties whose input has been used.
for v = the_parties{sets};% Now count the first parties too.
    usedinputs{1} = insert(v, usedinputs{1}, 0);
end
%If there is one set, then all of them are independent.
for n=1:N;
    if ismember(n,the_parties{1})&&ismember(n,the_parties{sets});
        fprintf('Party %d is causally independent', n);
        cprintf('blue', '  --- *** ---');
        fprintf('\n');
    end
end
term_input{N} = [];
if length(usedinputs{1})~= N;
    for u = find(~ismember([1:N],usedinputs{1}));% the parties that did not have an input
            if ismember(u,the_parties{1});%If that party is 'last' then this party is independent.
                fprintf('Party %d is causally independent', u);
                cprintf('blue', '  --- *** ---');
                fprintf('\n');
                % And it is better to consider this party as first, and also last (because markovianity looks at those to create terms).
                the_parties{sets} = insert(u, the_parties{sets},1);%this destroys the sorting so I'm fixing it later
            end
    end
end
the_parties{sets} = sort(the_parties{sets});
disp(['Time ' num2str(toc)]);

%***********************************************************************************************************
%**************************************** Markovianity *****************************************************
%***********************************************************************************************************

% For every party, find the parents and create a channel from all the outputs to the input. Then sort them.
% Terms and subterms are still useful,
% except for the arrows involving an input system connected to more than one output.
%  A = total_arrows(:,2);% the second column represents the parties of input systems involved in an arrow.
% Note that unique(A) not only displays the parties only once, but they are ordered too. 1 2 3 4... 
% [n, bin] = histc(A, unique(A));% n is the number of times an element of unique(A) is in A.
% Again, it refers to how many unique(A)(1) are in A, so also this info is ordered.
% All n shoud be one, except for those elements of unique(A) that appear multiple times in A.
% Bin, is the index 
% Therefore, if 4 appears two times, and this bin is specified in bin(3), then n(3) > 1.
% bin is the index of the bin in which each element of A belongs to.
% i.e. 1st element belongs to 3rd bin, 2nd element belongs to 2nd bin etc = [3 2 ...] 
% multiple = find(n > 1);% multiple finds the index of the elements of n: that denotes the number of elements appeared A within the range unique(A)(multiple)
% Again, bin, for each element of A indicates the index of the binrange = unique(A) it belongs to.
% Same elements in bin denotes same elements in A with the same index.\
% multiple denotes the index of an element in n that appears many times and is in the binrange at position multiple.

A  = total_arrows(:,2);% only the parties whose input has an arrow
tot = store;% input and output parties
bi = unique(A);
[n, bin] = histc(A, bi);
multiple = find(n > 1);
index    = find(ismember(bin, multiple));% The index of A's elements that appear multiple times.
sames = tot(index, :);% only the multiple input parties, with their outputs too.
manyarrows{N} = [];
a = 1;
sames1 = sames; % Buffer.



% Output system of the last parties
for u = the_parties{1};% The last parties.
    term{u,N+1} = traceoutcom(W, 2*u, dim);
end

multiarrows{N} = [];
multiarrows0{N} = []; % In the case of no subsystems, this never gets to be defined.
multiarrowsN{N} = [];
sys{N} = [];
subterm_tot{N,3} = [];% In dimension (:,1) is the term, in (:,2) the sys (whose complement is traced out 
% from W to obtain the term, in (:,3) is dim_mark that was used in traceoutcom(W, sys, dim_mark).
% Create the subterm.
% Here the difficulty is that we need the positions of all output systems or subsystems of (dim_mark) that are involved in the arrows towards partyin.
% Also, we need the position of party in. But with inserting subs, we loose track. So, in input all subs starting from the parties with the biggest number
% and like this subsequent parties subsystems or systems positions will not be altered.

part = bi(multiple);
if ~isempty(part);
    for st = 1:length(part);
for partyin = part(st);%Only those parties whose input is involved in multiple arrows.
    dim_mark = dim;
    multiarrows{partyin} = store(find(ismember(store(:,2),partyin)),:);% =[partyout partyin subs] with the same partyin
%     Now we create multiarrowsN{partyin} with the double appeared partyout partyin removed.
    C = multiarrows{partyin};%buffer
    c1 = C(:,1:2); %removed the info about the subsystems, to do the sorting
    [cuniq, cindex, cindexinuni] = unique(c1,'rows'); % cuniq has the unique sorted elements of c1, cindex has the indices of those elements in c1
%     cindex has the index of cuniq for every element of c1.
%     What we need is cindex to extract the unique (in terms of the first two elements of each row) elements from C = multiarrrows{paryin}
    multiarrowsN{partyin} = C(cindex,:);
    position = zeros(1,size(multiarrows{partyin},1));
acc = 0;
        for row1 = 1:size(multiarrowsN{partyin},1);% running the rows of subarrows 
            partyout = multiarrowsN{partyin}(row1,1);% The party whose output is linked with the input of party N.
            % Create dim_mark, dimensions subsystems inserted so that the right complement is traced out.
            i = i+1;
            if  multiarrowsN{partyin}(row1,3)~=0;% if the output party we're dealing with had subsystems, dim has to be altered 
                position(row1) = 2*partyout;% position in the dim
                dim_mark(position(row1)+acc) = []; % Remove the ouput dimensions to insert the subdim
                dim_mark_new{partyin,row1} = insert(subdim{partyout}, dim_mark, position(row1)+acc-1);
                dim_mark = dim_mark_new{partyin,row1};
% dim_mark keeps increasing. So the position of subsequent insert(subdim.. will have to shifted. 
% multiarrowsN is sorted so first the left subdims will be inserted, and keep going to the right. So the position increases by:
                acc = acc+length(subdim{partyout})-1;
            end
        end
        % Now all the subdims have been inserted succesfully 
        % Multiarrows has to include partyin and has to be sorted in ascend way now for the positions to be determined
    
    multiarrows0{partyin} = multiarrows{partyin};% to be used later on this form, because right now one line is added but is only needed in this loop
    [D1,D2] = sort(multiarrows0{partyin}(:,1));
    multiarrows0{partyin} = multiarrows0{partyin}(D2,:);
    
    multiarrows{partyin} = [multiarrows{partyin}; partyin partyin 0];%only useful for this loop
    [d1,d2] = sort(multiarrows{partyin}(:,1));
    multiarrows{partyin} = multiarrows{partyin}(d2,:);
    ls3 = 0;
     i = 0;
    for row = 1:size(multiarrows{partyin},1);% running the rows of subarrows 
            partyout = multiarrows{partyin}(row,1);% The party whose output is linked with the input of party N.
            % Find the position of this output (sub)system of this party
            i = i+1;
            if partyout ~= partyin;
                if ~isempty(subdim{partyout});
                    sys{partyin}(i) = 2*partyout - 1 + multiarrows{partyin}(row,3) + ls3;
                     ls3 = ls3 + size(subdim{partyout},2)-1;
                    if row>1 && partyout == multiarrows{partyin}(row-1,1);% if the partyout is the same as the previous then cancel the counting of ls3
                     sys{partyin}(i) = sys{partyin}(i) - size(subdim{partyout},2)+1 ;
                     ls3 = ls3 - size(subdim{partyout},2) +1;% because the ouput systems of the partyout==partyin are irrelevant, hence they should not be counted
                       
                    end
                else
                    sys{partyin}(i) = 2*partyout + ls3;
                end
            else
                sys{partyin}(i) = 2*partyin-1 + ls3;
                partysave = 2*partyin-1 +ls3;% this is the index in sys{partyin} of the input of partyin.
            end
    end
     
        if length(sys{partyin}) ~= size(multiarrows{partyin},1);
% Sanity check: # of systems whose complementary to be traced out == length of subarrows{partyin}
               cprintf('red', 'Something is wrong: traceoutcom has wrong input: sys');
        end
        
% Here all the output subsystems have been entered into dim_mark
% and if there weren't any, dim_mark is dim. So we're ready to create the subterm
        buff3 = traceoutcomN(W, sys{partyin}, dim_mark);% Unnormalized term
        norm3 = prod(dim_mark(sys{partyin}))/dim(2*partyin-1);% product of input and output systems involved/dim of the input
        subterm_tot{partyin,1} = buff3/trace(buff3)*norm3;% Normalized term.
        subterm_tot{partyin,2} = sys{partyin};% to keep track of the subsystems+systems involved (knowing the parties linked to partyin)
        subterm_tot{partyin,3} = dim_mark;
% subterm is a tensor of many outputs and one input. Like Ao Do Bo and Ci if ADB are sending to C. 
% The term would be a tensor of AoBoCiDo so we need to permute the tensor product such that the input system is last.
        if subterm_tot{partyin,2}(length(sys{partyin}))~=partysave;% if the last entry is not the partyin
%             dim_temp = zeros(1,length(sys{partyin}));
            spotin = find(sys{partyin}==partysave); % the input system has to go last in subeterm_tot{partyin,1}
            permusys = [1:spotin-1 spotin+1:length(sys{partyin}) spotin];
% dim_mark contains all the subsystems of parties involved in this multiarrow and sys{partyin} the indices of the relevant systems
            permudim = dim_mark(sys{partyin});
            subterm_tot{partyin,1} = syspermute(subterm_tot{partyin,1}, permusys, permudim);
            subterm_tot{partyin,2} = sys{partyin}(permusys);
            subterm_tot{partyin,3} = permudim;
        end
            
end
    end
end

% Some term{out,in} (those that has no output subsystems but go to an input subsystem) 
% has to be emptied (because subterm or subterm_tot{} is going to be used instead)

for line = 1:size(sames,1);
    if sames(line,3)==0;% if there is an arrow from sames(line,1) to sames(line,2) with no subsystems involved
    term{sames(line,1), sames(line,2)}=[];
    end % the other case is taken care of below        
end

all_arrows = [primal_arrows; secondary_arrows];
% Input system of first parties
for first = the_parties{sets};
        term_input{first} = traceoutcom(W, 2*first-1, dim);
end
% Input system of other parties that have no arrows coming to them
for no_input = 1:N;
    if ~ismember(no_input, all_arrows(:,2));
        term_input{no_input} = traceoutcom(W, 2*no_input-1, dim);
    end
end


% if max_sub==0;% If there are zero subsystems, subterm is not yet created. So create one
%     subterm{N,N,max_sub}=[];
% end

manyarrows = multiarrows0;
if n_subs>0;
%***************************************** Removing from subterms the manyarrows *****************************************
% multiarrows and subarrows have the modified labelling. However, subterm has the original one.
% Therefore, I need to change it to the right one whenever necessary.
for pr = 1:length(multiarrows);
    if ~isempty(manyarrows{pr});% there's multiple outputs linked to the input of pr.
              
            for gr = 1:size(manyarrows{pr},1);% the lines of subarrows
                partyo = manyarrows{pr}(gr,1);% the partyout
                if manyarrows{pr}(gr,3)~=0;% then there a subsystem involved
                    ai = manyarrows{pr}(gr,3);% the subsystem
                        if ~isempty(indexing{partyo})% if that party had an open subsystem
%                     then the ai subsystem we're checking is actually the find(indexing{partyo}==1)(ai)
                        ais = find(indexing{partyo}==1);
                            if ai<=length(ais);
                                ai = ais(ai);
                            end
                        end 
                        if ~isempty(subterm{manyarrows{pr}(gr,1), manyarrows{pr}(gr,2),ai});
                    subterm{manyarrows{pr}(gr,1), manyarrows{pr}(gr,2),ai} = [];
                        end
                end
            end
            
    end
end
end
                   
% The following are ready to construct W:
% The term{n2,n1} (from n2 to n1), term{n, N+1} which is identity on the last party n,
% term_input{n} is the input system of the first party n,
% subterm{n2,n1,p} (from subsystem p of n2, to n1, where n1 ~= any(sames(:,2))
% and subterm_tot{partyin,1} which is a channel from parties
% store(store(:,length(subarrows{partyin}))==partyin, 1) to partyin

%******************************************* Constructing Wtest ********************************************


Wtest = 1;
% order is a row keeping track of the order of the systems [ 1 2 1 4 2 3 3 4 ]
order{1} = []; % 2*N (# systems) + tot_subs (# total subsystems) - n_subs (# output systems that have subsystems)
order{2} = []; % The dimensions of each input/output system that enters
order{3} = []; % 0 if it's an input system, 1 if it's an output system for each element of order{1} of dimensions order{2}.
order{4} = []; % Here the labelling of the subsystem used will be stored. Taken from subarrows{}(:,3).
for pp = 1:N;
     % if no single links are going to pp, it might be because several links do
            if ~isempty(subterm_tot{pp,1});% if there is a subarrow to party pp
                Wtest = tensor(Wtest, subterm_tot{pp,1});
                for linesub = 1:size(manyarrows{pp},1);% running for the parties that send to party pp
                    pout = manyarrows{pp}(linesub,1);
                    order{1} = [order{1} pout];
                    order{3} = [order{3} 1];
                    if ~isempty(subdim{pout})
                        order{2} = [order{2} subdim{pout}(manyarrows{pp}(linesub,3))]; % dim of the correct subsystem that goes to pp

                        if ~isempty(indexing{manyarrows{pp}(linesub,1)}); % the right labelling of the subs needs to written in order{4}
                            label = find(indexing{manyarrows{pp}(linesub,1)}==1);
                            order{4} = [order{4} label(manyarrows{pp}(linesub,3))];
                        else
                            order{4} = [order{4} manyarrows{pp}(linesub,3)];
                        end
                        
                    else
                        order{2} = [order{2} dim(2*pout)];% dim if the output system of pout.
                        order{4} = [order{4} 0];
                    end
                end
                order{1} = [order{1} pp];
                order{2} = [order{2} dim(2*pp-1)]; % dim of the input system of pp.
                order{3} = [order{3} 0];
                order{4} = [order{4} 0];
%                     fprintf('three');
%                     fprintf('\n')
            end
            % if pp is last or if is first
            
            if ~isempty(term_input{pp});% if it's first
                Wtest = tensor(Wtest,term_input{pp});
                order{1} = [order{1} pp];
                order{2} = [order{2} dim(2*pp-1)];
                order{3} = [order{3} 0];
                order{4} = [order{4} 0];
%                 fprintf('zero');
%                 fprintf('\n')
            end
            if ~isempty(term{pp,N+1});% if it's last
                Wtest = tensor(Wtest,term{pp, N+1});
                order{1} = [order{1} pp];
                order{2} = [order{2} dim(2*pp)];
                order{3} = [order{3} 1];
                order{4} = [order{4} 0];
%                 fprintf('four');
%                 fprintf('\n')
            end
            
            
    for qq = 1:N;
        if qq ~= pp;
            if ~isempty(term{pp,qq});%if the output is connected to an input
                Wtest = tensor(Wtest, term{pp, qq});
                order{1} = [order{1} pp qq];
                order{2} = [order{2} dim(2*pp) dim(2*qq-1)];% dims of output of qq, input of qq
                order{3} = [order{3} 1 0];
                order{4} = [order{4} 0 0];
%                 fprintf('two');% helps debugging
%                 fprintf('\n')
            else
                if max_sub;
                if size(subterm,3)~=1;
                for mi = 1:size(subterm,3);
                    if ~isempty(subterm{pp,qq,mi});
                        if ~isempty(indexing{pp});
                            if indexing{pp}(mi)==1;% this means the system was not open
                        order{2} = [order{2} subdimor{pp}(mi) dim(2*qq-1)];
                        Wtest = tensor(Wtest, subterm{pp,qq,mi});
                        order{1} = [order{1} pp qq];
                        order{3} = [order{3} 1 0];
                        order{4} = [order{4} mi 0];
%                         fprintf('x_one');% helps debugging
%                         fprintf('\n')
                            end
                        else
                        Wtest = tensor(Wtest, subterm{pp,qq,mi});
                        order{1} = [order{1} pp qq];
                        order{2} = [order{2} subdimor{pp}(mi) dim(2*qq-1)];% dims of ouput of the subsys mi of pp, input of qq
                        order{3} = [order{3} 1 0];
                        order{4} = [order{4} mi 0];
%                         fprintf('one');% helps debugging
%                         fprintf('\n')
                        end
                    end
                end
                end
                end
            end
        end
    end
end
% At this point the Wtest has all the terms it should have. 
% Therefore, in 'order' there are the pairs in the order they appear on the constructed Wtest.

orderline1 = order{1};
orderline3 = order{3};
orderline4 = order{4};

% First of all, if manyarrows is nonempty, then reshuffling is needed before checking anything else.
% To do that, a dim_final has to be created that has all the dimension of all subsystems.
dim_final{1} = dim;
ls4 = 0;
if ~cellfun(@isempty, subdim);
    for mm = find(~cellfun(@isempty, subdim)); % the index of non-empty subdim / the party whose output has subsystems.
        dim_final{1}(2*mm + ls4) = [];
        dim_final{1} = insert(subdim{mm}, dim_final{1}, 2*mm+ls4-1);
        ls4 = ls4+length(subdim{mm})-1;% That's how many extra elements now dim_final has.
    end
end
dim_f = dim_final{1};
keep{1} = [];
dims = order{2};
savett = [];
savetj = [];
lend = length(dims);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reshuffling for multi arrows. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now dim_f has the dimensions of all output subsystems. Time to reshuffle things. First, deal with the multiple arrows. It may be enough. 
% If there are multiple arrows to one party.
for ww = find(~cellfun(@isempty, manyarrows)); % the index of non-empty subarrows == the party that has multiple arrows towards them.
            
for tj = size(manyarrows{ww},1):-1:1;%the position needs to change from ww-3 to ww-2 to ww-1 for the parties that output to party ww
    inout = find(order{1}==ww); num = find(order{3}(inout)==0);%there is only one 0 (input) at all the places of inout)
%     so num is the position of the input of ww no matter how many outputs have been before that in order{1}
     savetj = [savetj tj];

              tt = inout(num)-tj;% positions: from the first to the last party that outputs to party ww
                savett = [savett tt];
%        Sanity check, there must be an output at this position 
                if order{3}(tt);% 1 as there should be an output
%        There must be consecutive ones corresponding to the multi-arrows.
%        Time to reshuffle: 
%        order{1}(tt) needs to go next to the first element of find(order{1} == order{1}(tt))
                    p_out = order{1}(tt);% the party that send out on party ww
%                     Here we're looking for the position of the input of party p_out so that we place it next to it.
%                     But because we don't know if the output of p_out is before or after its input,
%                     we have to look for its input which will be where order{1}==p_out but also where order{3}==0
    
                    F = find(order{1}==p_out);
                    d1 = F(find(order{3}(F)==0)); % d1 is the position of the input of party p_out.                   
                    leng = 1;
                    
                    if ~isempty(subdimor{p_out});
                        leng = length(subdimor{pout});
                    end
                        
                    if tt > d1+leng; % then the output needs to go to the LEFT
                        
%                     fprintf('g - ');% helps debugging
%                     fprintf('\n')
                    savedim = [prod(dims(1:d1)) prod(dims(d1+1:tt-1)) dims(tt) prod(dims(tt+1:length(dims)))];
                    Wtest = syspermute(Wtest, [1 3 2 4], savedim);
                    
                    orderkeep{p_out} = order{1};
                        for yy = 1:4;
                        order{yy} = order{yy}([1:d1, tt,d1+1:tt-1, tt+1:lend]);
                        end
                    dims = order{2};
                    
                    elseif tt < d1  % then the output needs to go to the RIGHT 
                     
%                     fprintf('g - reverse'); % helps debugging
%                     fprintf('\n')
                    dimdim = [prod(dims(1:tt-1)) dims(tt) prod(dims(tt+1:d1)) prod(dims(d1+1:lend))];
                    Wtest = syspermute(Wtest, [1 3 2 4], dimdim);
                    
                    orderkeep{p_out} = order{1};
                        for yy = 1:4;
                        order{yy} = order{yy}([1:tt-1, tt+1:d1, tt, d1+1:lend]);
                        end
                      dims = order{2};
                   end
                end

end             
end
ordertest{4} = [];
ordertest{1} = [];
%***************************************** Builiding ordetest{1} *****************************************
for ii = 1:N;
    ordertest{1} = [ ordertest{1} ii ii ];% recording the input and output system of party ii
    if ~isempty(subdim{ii})&& ~ismember(ii,the_parties{1});% if there are subs and the party is not last
        for oo = 1:length(subdim{ii})-1;% -1 because one output system has already been recorded two lines above
        ordertest{1} = [ordertest{1} ii];
        end
    end
end

%***************************************** Reshuffling for regular arrows. *****************************************
% Here the subdim{} that is not empty matters. However, there might be subsystems that were not counted as open because the
% parties that had them were last. So I'll empty the subdim{} of the parties that were last : the_parties{1}

for lastp = the_parties{1};
subdim{lastp} = [];
end

% Info about those is stored in store where the last column is 0
% So this section is finding the input, then its (sub)ouput and is placing it next to the input
% Aim: order{1} = 1 1 2 2 3 3 3 4 4 etc.]
store(find(store(:,1)==0),:,:) = []; % Remove the zero rows.
arrows = store(find(store(:,3)==0),:); % Now arrows have the simple arrows info
wtf = order{1};
dimss = order{2};
if sum(order{1}==ordertest{1}) ~= length(order{1});% if order{1}~=ordertest{1}
for parin = 1:N;
    
    spotin = find(order{1}==parin);% spots of the input and output systems of parin
    findinput = find(order{3}(spotin)==0); % the position of the input in the above spots
    spotin = spotin(findinput);% Take the spot of input in order{1}
    
    if isempty(subdim{parin});
        spotout = find(order{1}==parin);% same story
        findoutput = find(order{3}(spotout)==1); % same
        spotout = spotout(findoutput); % Spot of the output system on order{1}
    
        targetout = spotin+1;% target of the output system
            spotarget{parin} = [spotout targetout];
        if targetout<spotout;
            dimt2s = [prod(dimss(1:targetout-1)) prod(dimss(targetout:spotout-1)) dimss(spotout) prod(dimss(spotout+1:length(dimss)))];
            Wtest = syspermute(Wtest, [1 3 2 4], dimt2s);
        
            for yy = 1:4;
                order{yy} = order{yy}([1:targetout-1, spotout, targetout:spotout-1, spotout+1:length(dimss)]);
            end
            dimss = order{2};
            
        elseif spotout<targetout; % if spotout==targetout nothing will happen
            dims2t = [prod(dimss(1:spotout-1)) dimss(spotout) prod(dimss(spotout+1:targetout-1)) prod(dimss(targetout:length(dimss)))];
            Wtest = syspermute(Wtest, [1 3 2 4], dims2t);
        
            for yy = 1:4;
                order{yy} = order{yy}([1:spotout-1, spotout+1:targetout-1, spotout, targetout:length(dimss)]);
            end
            dimss = order{2};
            
        end
        
    else
% if this party has subsystems that need to be put next to the input
% But every time a subsystem goes into the right position,
% the position of the rest of the subsystems and the input system changes.
% This will be the case every time spotout<targetout. So I shift 
    
        n_Subs = length(subdim{parin});
        spotout = find(order{1}==parin);% spot input and outputs 
        findoutput = find(order{3}(spotout)==1); % indices of the outputs 
        spotouts = spotout(findoutput);%
        for jk = 1:n_Subs;% for every output system
        
            spotout = spotouts(jk);
            targetout = spotin+1;
%             Same story as above
            if targetout<spotout;
            dimt2s = [prod(dimss(1:targetout-1)) prod(dimss(targetout:spotout-1)) dimss(spotout) prod(dimss(spotout+1:length(dimss)))];
            Wtest = syspermute(Wtest, [1 3 2 4], dimt2s);
        
                for yy = 1:4;
                order{yy} = order{yy}([1:targetout-1, spotout, targetout:spotout-1, spotout+1:length(dimss)]);
                end
            dimss = order{2};
            
            elseif spotout<targetout; % if spotout==targetout nothing will happen
               
            dims2t = [prod(dimss(1:spotout-1)) dimss(spotout) prod(dimss(spotout+1:targetout-1)) prod(dimss(targetout:length(dimss)))];
            Wtest = syspermute(Wtest, [1 3 2 4], dims2t);
        
                for yy = 1:4;
                order{yy} = order{yy}([1:spotout-1, spotout+1:targetout-1, spotout, targetout:length(dimss)]);
                end
                dimss = order{2};
                spotouts(spotouts<spotin)=spotouts(spotouts<spotin)-1;
                spotin = spotin-1;
                
                
            end
          
        end
    end
            
            
end
end

% The above guarantees that order{1}= [1 1 1 2 2 5 5 4 4 3 3] so the outputs will be next to the inputs but 
% the order of the parties might still not be the right one. So first I create ordertest{1} and then compare it to order{1} 
% and re-shuffle order{1} when necessary

 fprintf('\n')    
for yy = 1:4;
    order{yy};
end
%  Now check if order{1} = [1 1.. 2 2.. 3 3 ... etc]
%  Create the matrix, to compare it.
%  In these ordertest{i} matrices, I do not count the output subsystems of the last parties.
%  Because they were not counted in the building of the order{i}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reshuffling of whole parties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subdimspe = subdim;
for lastparties = the_parties{1};
subdimspe{lastparties} = [];
end

count2 = 0;
singlesub = zeros(1,N);
dims = order{2};

if sum(order{1}==ordertest{1}) ~= length(order{1});
for partyno = 1:N;
    if ~isempty(subdimspe{partyno});
        len = length(subdimspe{partyno})-1;
        % quick fix for later
        if len==0;
            singlesub(partyno) = 1; % if the subsystems have only one system left (because the rest were open)
        end
    else
        len = 0;
    end
    
    order1 = order{1}(2*partyno-1+count2:2*partyno+count2+len);% this should be 111 or 2222 or 33 etc i.e. [input outputs]of the same party
    ordertest1 = ordertest{1}(2*partyno-1+count2:2*partyno+count2+len);%this _is_ the above explanation.
    count2 = count2+len;
    if sum(order1==ordertest1)~=length(order1); % if they are not the same
% Reshuffling is needed
% ordertest{1}(count):ordertest{1}(count+1+len) has the right numbers in
% Need to find where those numbers are in order{1} (they have the right order "input outputs" at least)
% and put them in the right position within order{1}
        spot = find(order{1}==partyno);% position of input and outputs
        ind = find(order{3}(spot)==0); % position of input in the above
        spot = spot(ind); % the position in order{1} where the input of the partyno is
     
        target = find(ordertest{1}==partyno,1); % this is where the order{1}(spot:spot+1+len) has to go
if target<spot % sometimes it won't be the case because the subsystems have not been sorted yet.
        prod1 = prod(dims(1:target-1));
   
        prodn = prod(dims(spot+1+len+1:length(dims)));
    
    dimswap = [prod1 prod(dims(target:spot-1)) prod(dims(spot:spot+1+len)) prodn];
    Wtest = syspermute(Wtest, [1 3 2 4], dimswap);
    swaps{partyno} = [spot target dimswap order{1} ];
    
%      target<spot so the element has to go left 
%     Now all the orders need to be swapped too

    for yy = 1:4;
        order{yy} = order{yy}([1:target-1, spot:spot+1+len, target:spot-1, spot+1+len+1:length(dims)]);
    end
dims = order{2} ;

end
    end
end 
end
       
     
if sum(ordertest{1} == order{1}) == length(ordertest{1});
    % Now we need to test the order of the subsystems
    for aa = 1:N;
        ordertest{4} = [ordertest{4} 0]; % for the input system
        
        if ~isempty(subdim{aa})&& ~ismember(aa,the_parties{1});% if there are subs and aa is not last
            
                if ~isempty(indexing{aa});
                    for dd = find(indexing{aa}==1);% the names of the subs that were not open 
                        ordertest{4} = [ordertest{4} dd];
                    end
                else % If there were no open output subs in aa, then take the labelling from subdim{aa}
                for ss = 1:length(subdim{aa});
                    ordertest{4} = [ordertest{4} ss];                    
                end
                end
        
        else
            ordertest{4} = [ordertest{4} 0]; % for the output system
        end
    end
end
    
% Now compare 
    
% Info about subdims is in subdimspe now so update the following
subs = sum(cellfun(@sum, subdimspe));% Is the sum of all the dimensions of subsystems
max_sub = max(cellfun(@length, subdimspe));% the number of max subsystems composing one system
tot_subs = sum(cellfun(@length, subdimspe)); % the number of all subsystems
n_subs = length(find(~cellfun(@isempty, subdimspe))); % the number of output systems (parties) that are divided to subsystems
% Now singlesub contains zeros except on parties that have one subsystem only left, in which case it's one.
% In that case, order{4} does not need any reshuffling of subsystems so I can put subdimspe{singlesub} = [];

for V = find(singlesub);
subdimspe{V} = [];
end

if max_sub>1;

jjj = 0;
if sum(order{4} == ordertest{4})~= length(order{4}); % Reshuffling is needed.
    jjj = 1;% checking

subs = sum(cellfun(@sum, subdim));% Is the sum of all the dimensions of subsystems
max_sub = max(cellfun(@length, subdim));% the number of max subsystems composing one system
tot_subs = sum(cellfun(@length, subdim)); % the number of all subsystems
n_subs = length(find(~cellfun(@isempty, subdim))); % the number of output systems (parties) that are divided to subsystems

    for n_psub = 1:n_subs; % the number of parties that have ouput subs and might need reshuffling
        for psub = find(~cellfun(@isempty,subdim));% each party 
%             The outputs we're looking for are sitting right next to the input of the parties. They only need rearrangement.
            psub_in = find(order{1}==psub,1);% position of the input
            psub_out = psub_in+1;
            if sum(order{4}(psub_out:psub_out+length(subdim{psub})-1)==ordertest{4}(psub_out:psub_out+length(subdim{psub})-1))~=length(subdim{psub});
                take = order{4}(psub_out:psub_out+length(subdim{psub})-1);
                [sorted, permutati] = sort(take);
%  For the real permutation I need 1 2 3 ..until psub_in then permutati + psub_in then psub_in+length(permutati): end

                permutation = [1:psub_in permutati+psub_in psub_in+length(permutati)+1:length(order{2})];
                Wtest = syspermute(Wtest, permutation, order{2});
%                 fprintf('swapping of subs happened'); % helps debugging
%                 fprintf('\n')
                for yy = 1:4;
                    order{yy} = order{yy}(permutation);
                end
            end
        end
    end
 
end
end

%  Now check 

if order{4} == ordertest{4};
    diff = abs(Wtest - W);
        if isempty(find(diff>epsilon,1));
            cprintf('blue', 'The process is Markovian');
              fprintf('\n')
        else
            cprintf('cyan', 'The process is _not_ Markovian');
              fprintf('\n')
            diff_elements = length(find(diff>epsilon))% these is the number of the non-zero elements of the matrix abs(Wtest-W)
        end
end

for yy = 1:4;
    order{yy};
end

% Output the DAG
R = all_arrows;
[f, g, h] = unique(R,'rows');% Now f has all the unique arrows of R
DG = sparse(f(:,1), f(:,2), true, N,N);
view(biograph(DG))
