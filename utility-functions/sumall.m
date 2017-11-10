function a = sumall(a, n)
% FORMAT s = sumall(a)
% Sums all elements of a.
% Very efficient in terms of speed and memory for a number of non 
% singelton dimensions <= 10. 
% Try to be not too ineficient for higher dimensions
    if isempty(a)
        a = 0;
        return
    end
    if nargin < 2
        n = sum(size(a) > 1);
    end
    switch n
        case {0, 1}
            a = sum(squeeze(a));
        case 2
            a = sum(sum(squeeze(a)));
        case 3
            a = sum(sum(sum(squeeze(a))));
        case 4
            a = sum(sum(sum(sum(squeeze(a)))));
        case 5
            a = sum(sum(sum(sum(sum(squeeze(a))))));
        case 6
            a = sum(sum(sum(sum(sum(sum(squeeze(a)))))));
        case 7
            a = sum(sum(sum(sum(sum(sum(sum(squeeze(a))))))));
        case 8
            a = sum(sum(sum(sum(sum(sum(sum(sum(squeeze(a)))))))));
        case 9
            a = sum(sum(sum(sum(sum(sum(sum(sum(sum(squeeze(a))))))))));
        case 10
            a = sum(sum(sum(sum(sum(sum(sum(sum(sum(sum(squeeze(a)))))))))));
        otherwise
            a = sum(a(:));
    end
end
