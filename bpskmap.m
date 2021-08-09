
function [y,A] = bpskmap(x)
            
            A = [-1,1]; %Constellation symbol Alphabet

            y = A(1)*(x==1) + A(2)*(x==0); % symbol vector
end