function HD = HDiff(A,B)
    sA = length(A);
    sB = length(B);
    daB = zeros(sA,1);
    dbA = zeros(sB,1);
    for i = 1:sA
        daB(i) = min(abs(B-A(i)));
    end
    for i = 1:sB
        dbA(i) = min(abs(A-B(i)));
    end
    HD = max([daB;dbA]);
end