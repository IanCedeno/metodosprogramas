function [x, iter] = jacobi_matricial_parte2(A, b, tolerancia, max_iter)
    n = length(b);
    x = zeros(n, 1);
    iter = 0;

    fprintf('Iteración    x(k)         y(k)         z(k)         Er,x         Er,y         Er,z\n');

    while iter < max_iter
        x_prev = x;
        for i = 1:n
            sum = 0;
            for j = 1:n
                if j ~= i
                    sum = sum + A(i, j) * x_prev(j);
                end
            end
            x(i) = (b(i) - sum) / A(i, i);
        end

        % Calcular errores relativos porcentuales
        errors = 100 * abs(x - x_prev) ./ abs(x);

        fprintf('%d            %.8f    %.8f    %.8f    %.6f%%    %.6f%%    %.6f%%\n', iter, x(1), x(2), x(3), errors(1), errors(2), errors(3));

        iter = iter + 1;

        if max(errors) < tolerancia
            break;
        end
    end
end
