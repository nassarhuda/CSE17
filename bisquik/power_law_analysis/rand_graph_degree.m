function A = rand_graph_degree(deg,ntrials)
% RAND_GRAPH_DEGREE Generate a random sample of a graph with a prescribed
% degree seqeunce or a perturbed version.

mypath = fileparts(mfilename('fullpath'));
deg = deg(:); % ensure the column dimension

write_degree_file('pltest.degs',deg);

trial = 1;
cmd = sprintf('%s pltest.degs -f pltest', fullfile(mypath,'..','bisquik-master','bisquik'));
status = system(cmd);
if status==255
    fprintf('Try %4i - Not a graphical sequence\n',trial);
    trial = trial + 1;
    while trial <= ntrials
        deg2 = deg + randi([0,1],numel(deg),1);
        if mod(sum(deg2),2) ~= 0
            deg2(end) = deg2(end)+1;
        end
        write_degree_file('pltest.degs',deg2);
        
        status = system(cmd);
        if status ~= 255
            break;
        end
        fprintf('Try %4i - Not a graphical sequence\n',trial);
        trial = trial + 1;
    end
    
    if status == 127
        % we never got a successful sample
        % TODO improve this error case
        assert(false);
    end
end

A = read_edges('pltest.edges');