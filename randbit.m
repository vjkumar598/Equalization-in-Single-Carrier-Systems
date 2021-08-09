%Random Bit Generator
%Bits are generated with equal probability,i.e uniformly distributed

function y = randbit(n)
            y = round(rand(1,n));
end
    