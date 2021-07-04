function idc = whatever()
    hello = zeros(6,1);
    myNest
    function bird = myNest()
        bird = 1;
        hello(2) = 100;
    end
    disp(hello)
    idc = 4;
end