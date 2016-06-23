
function a = forward_backward_substitution(U,Y)
    a = forward(U.',Y);
    a = backward(U,a);
end

function a = forward(UT,Y)
    n = length(Y);
    a = zeros(n,1);
    for i = 1:n
        v = Y(i);
        for j = 1:i-1
            v = v-UT(j,i)*a(j);
        end
        a(i) = v/UT(i,i);
    end
end

function a = backward(U,a)
    n = length(a);
    for i = n:-1:1
        v = a(i);
        for j = n:-1:i+1
            v = v-U(i,j)*a(j);
        end
        a(i) = v/U(i,i);
    end
end
