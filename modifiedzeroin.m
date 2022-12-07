function [root, info] = modifiedzeroin3036082191(func, Int, params)
    func_tol = params.func_tol; root_tol = params.root_tol; maxit = params.maxit;
    a = Int.a; b = Int.b; f = func;
    
    iter = 0;
    interpIter = 0;
    
    x0 = a; x1 = b; x2 = (a + b) / 2;
    
    while iter <= maxit
        f0 = f(x0); f1 = f(x1); f2 = f(x2);
        if x1 - x0 < root_tol
            if abs(f2) < func_tol
                root = x2;
                info = 0;
            else
                root = x2;
                info = 1;
            end
            return
        end
        
        y0 = f0; y1 = f1; y2 = f2;
        x3 = (y1 * y2 * x0) / ((y0 - y1) * (y0 - y2)) + (y0 * y2 * x1) / ((y1 - y0) * (y1 - y2)) + (y0 * y1 * x2) / ((y2 - y0) * (y2 - y1));
        f3 = f(x3);
        interpIter = interpIter + 1;
        
        try
            if abs(f3) < func_tol
                root = x3;
                info = 0;
                return
            end
        catch
        end
        
        if (x3 < x0) || (x3 > x1)
            fl = y0;
            diff = (x1 - x0) / 2;
            fm = f2;
            bisectionIter = 1;
            while (abs(fm) > func_tol) && (diff > root_tol) && (bisectionIter < 12)
                bisectionIter = bisectionIter + 1;
                if fl * fm < 0
                    x1 = x2;
                else
                    x0 = x2;
                    fl = fm;
                end
                x2 = (x1 + x0) / 2;
                diff = (x1 - x0) / 2;
                fm = f(x2);
            end
            interpCheck = fm; interpIter = 0;
        else
            if x3 < x2
                if f3 * f2 < 0
                    x0 = x3; x1 = x2;
                    fl = f3;
                elseif f3 * f0 < 0
                    x1 = x3;
                else
                    x0 = x2;
                    fl = f2;
                end
            else
                if f3 * f2 < 0
                    x0 = x2; x1 = x3;
                    fl = f2;
                elseif f3 * f1 < 0
                    x0 = x3;
                    fl = f3;
                else
                    x1 = x2;
                end
            end
            
            x2 = (x0 + x1) / 2;
            
            if interpIter == 1
                interpCheck = f3;
            end
        end
        
        if (interpIter > 7) && (abs(f3) > abs(interpCheck / 2))
            diff = (x1 - x0) / 2;
            fm = f(x2);
            bisectionIter = 1;
            while (abs(fm) > func_tol) && (diff > root_tol) && (bisectionIter < 10)
                bisectionIter = bisectionIter + 1;
                if fl * fm < 0
                    x1 = x2;
                else
                    x0 = x2;
                    fl = fm;
                end
                x2 = (x1 + x0) / 2;
                diff = (x1 - x0) / 2;
                fm = f(x2);
            end
            interpIter = 0;
        end
        
        iter = iter + 1;
    end
    
    if abs(f3) < func_tol
        root = x3;
        info = 0;
    else
        root = x3;
        info = 1;
    end
    return
end