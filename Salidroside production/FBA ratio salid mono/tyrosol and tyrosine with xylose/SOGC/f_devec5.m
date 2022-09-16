function [bestmem , bestval, nfeval , iter] = f_devec5 ( fname , VTR , D ,...
    XVmin , XVmax, y , NP , itermax , F , CR , strategy , refresh ,...
    resume)
% -------------------------------------------------------------------------

%-----Initialize population and some arrays-------------------------------
pop = zeros(NP,D); %initialize pop to gain speed
%----pop is a matrix of size NPxD. It will be initialized-------------
%----with random values between the min and max values of the---------
%----parameters-------------------------------------------------------

if resume
    
    load('Last_DE_data','Last_DE_data');
    pop = Last_DE_data.pop;
    nfeval = Last_DE_data.nfeval;
    iter = Last_DE_data.iter+1;
    
else
    
    for i=1:NP
        pop(i,:) = XVmin + rand(1,D).*(XVmax - XVmin);
    end

    nfeval = 0; % number of function evaluations
    iter = 1;
    
end

pop=round(100*pop)/100;
% pop(:,end)=round(pop(:,end));

val      =  zeros(1,NP);          % create and reset the "cost array"

%------Evaluate the best member after initialization----------------------

ibest    = 1;           % start with first population member
val(1)   = feval(fname,pop(ibest,:),y);
bestval  = val(1);                 % best objective function value so far
worstval = val(1);                 % worst objective function value so far
nfeval   = nfeval + 1;

for i=2:NP                        % check the remaining members
    val(i) = feval(fname,pop(i,:),y);
    nfeval  = nfeval + 1;
    if (val(i) < bestval)           % if member is better
        ibest   = i;                 % save its location
        bestval = val(i);
    elseif (val(i) > worstval)
        worstval = val(i);
    end
    
end

bestmemit = pop(ibest,:);           % best member of current iteration

bestmem = bestmemit;                % best member ever

%------DE-Minimization---------------------------------------------
%------popold is the population which has to compete. It is--------
%------static through one iteration. pop is the newly--------------
%------emerging population.----------------------------------------

bm   = zeros(NP,D);              % initialize bestmember matrix
ui   = zeros(NP,D);        % intermediate population of perturbed vectors
rot  = (0:1:NP-1);               % rotating index array (size NP)
rotd = (0:1:D-1);                % rotating index array (size D)

%------------------------------------
while ((iter < itermax) && (abs(bestval - worstval) > VTR))
    popold = pop;                   % save the old population
    
    ind = randperm(4);              % index pointer array
    
    a1  = randperm(NP);             % shuffle locations of vectors
    rt  = rem(rot+ind(1),NP);        % rotate indices by ind(1) positions
    a2  = a1(rt+1);                 % rotate vector locations
    rt  = rem(rot+ind(2),NP);
    a3  = a2(rt+1);
    rt  = rem(rot+ind(3),NP);
    a4  = a3(rt+1);
    rt  = rem(rot+ind(4),NP);
    a5  = a4(rt+1);
    
    pm1 = popold(a1,:);             % shuffled population 1
    pm2 = popold(a2,:);             % shuffled population 2
    pm3 = popold(a3,:);             % shuffled population 3
    pm4 = popold(a4,:);             % shuffled population 4
    pm5 = popold(a5,:);             % shuffled population 5
    
    for i=1:NP                   % population filled with the best member
        bm(i,:) = bestmemit;     % of the last iteration
    end
    
    mui = rand(NP,D) < CR;   % all random numbers < CR are 1, 0 otherwise
    
    if (strategy > 5)
        st = strategy-5;		  % binomial crossover
    else
        st = strategy;		  % exponential crossover
        mui=sort(mui');	          % transpose, collect 1's in each column
        for i=1:NP
            n=floor(rand*D);
            if n > 0
                rtd = rem(rotd+n,D);
                mui(:,i) = mui(rtd+1,i); %rotate column i by n
            end
        end
        mui = mui';			  % transpose back
    end
    mpo = mui < 0.5;                % inverse mask to mui
    
    if (st == 1)                      % DE/best/1
        ui = bm + F*(pm1 - pm2);        % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 2)                  % DE/rand/1
        ui = pm3 + F*(pm1 - pm2);       % differential variation
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 3)                  % DE/rand-to-best/1
        ui = popold + F*(bm-popold) + F*(pm1 - pm2);
        ui = popold.*mpo + ui.*mui;     % crossover
    elseif (st == 4)                  % DE/best/2
        ui = bm + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;           % crossover
    elseif (st == 5)                  % DE/rand/2
        ui = pm5 + F*(pm1 - pm2 + pm3 - pm4);  % differential variation
        ui = popold.*mpo + ui.*mui;            % crossover
    end
    % ----------------------------------
    %        for i=1:NP
    %          for j=1:D
    %              if ui(i,j)<XVmin(j)
    %                  ui(i,j)=XVmin(j);
    %              elseif ui(i,j)>XVmax(j)
    %                  ui(i,j)=XVmax(j);
    %              end
    %          end
    %        end
    % ---------------------------------
    for i=1:NP
        for j=1:D
            if ui(i,j)<XVmin(j)||ui(i,j)>XVmax(j)
                ui(i,j)=XVmin(j) + rand(1)*(XVmax(j) - XVmin(j));
            end
        end
    end
    % ---------------------------------
    ui=round(100*ui)/100;
%     ui(:,end)=round(ui(:,end));
    % ---------------------------------
    %-----Select which vectors are allowed to enter the new population----
    for i=1:NP
        tempval = feval(fname,ui(i,:),y);   % check cost of competitor
        nfeval  = nfeval + 1;
        if (tempval <= val(i))  % if competitor is better than
            %                                 value in "cost array"
            pop(i,:) = ui(i,:);  % replace old vector with new one
            %                                (for new iteration)            
            val(i)   = tempval;  % save value in "cost array"
            
            %----we update bestval only in case of success to save time---
            if (tempval < bestval)     % if competitor better than
                %                                        the best one ever
                bestval = tempval;      % new best value
                bestmem = ui(i,:);      % new best parameter vector ever
            elseif (tempval > worstval)
                worstval = tempval;      % new worst value
            end
        end
    end %---end for imember=1:NP
    
    bestmemit = bestmem;       % freeze the best member of this
    %                                iteration for the coming
    % iteration. This is needed for some of the strategies.
    
    %----Output section---------------------------------------------------
    
    if (refresh > 0)
        if (rem(iter,refresh) == 0)
            fprintf(1,...
                'Iteration: %d,  Best: %f,  F: %f,  CR: %f,  NP: %d\n',...
                iter,bestval,F,CR,NP);
            for n=1:D
                fprintf(1,'best(%d) = %f\n', n , bestmem(n));
            end
        end
    end
    
    %------------------------
    Last_DE_data.pop = pop;
    Last_DE_data.nfeval = nfeval;
    Last_DE_data.iter = iter;
    save('Last_DE_data','Last_DE_data');
    %------------------------
    iter = iter + 1;
    %------------------------
end %---end while ((iter < itermax) ...
end