function esc=Esc()
% Set escape character separating folders
machine=computer;
if strcmp(machine(1:3),'MAC') 
    esc='/';
elseif strcmp(machine(1:3),'GLN') % unix cluster GLNXA64
    esc='/';
else % windows
    esc='\';
end

