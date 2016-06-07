function [W,X]= TransformLearning(W,ZCA,Y,iter,l2,l3,l4,p,step,cg_iter,thr,T0,debug,visual,saved_path)

    [M,n]   = size(W);
    X       = zeros(M,size(Y,2));
    YYT     = Y*Y';

    for i=1:iter
        if T0 == 0
            X   = soft(W*Y,thr);  
        else
            tmp     = W*Y;
            sorted  = sort(abs(tmp),'descend');
            X       = tmp.*bsxfun(@ge,abs(tmp),repmat(sorted(T0,:),M,1));
        end
        
        XYT = X*Y';
        for j = 1:cg_iter
            ZZ  = (W*W').^(p-1);
            cc  = -2*XYT + 2*W*YYT';                % 2WYY'-2XY'
            if M>n                                  % grad of log det
                cc2     =   (-2)*W*(inv(W'*W));     % over complete  
            else
                cc2     =   (-2)*inv(W*W')*W;       % less complete
            end
            cc3     =   2*W;                        %grad of ||W||_F
            cc4     =   2*p*((ZZ*W) - (((diag(ZZ))*(ones(1,n))).*W));   % grad of incoherence
            deD     =   l2*cc2 + cc + l3*cc3  + l4*cc4;
        
            if(j==1)
                g   =   deD;
                d   =   -g;
            else
                be  =   ((norm(deD,'fro'))^2)/((norm(g,'fro'))^2);
                g   =   deD;
                d   =   - g + be*d;
            end
            
            if debug == 1
                if M>n
                    logdet = -log(det(W'*W));
                else
                    logdet = -log(det(W*W'));
                end
                norm_g = sprintf('Norm of grad:%f log(det):%f',norm(g),logdet);
                disp(norm_g)
            end

            W = W + (step)*d;
        end
        
        %post-normalization of the rows of the transform
        SL          = diag(1./sqrt(sum((W').^2)));
        W           = SL*W;
        
        err_str     = sprintf('%d-th iteration,Err of reconstruction:%f',i,norm(W*Y-X,'fro'));
        disp(err_str);
        
        if visual == 1
            visualization     = displayDictionaryElementsAsImage((W*ZCA)', sqrt(M), sqrt(M),sqrt(n),sqrt(n));
            imwrite(visualization,[saved_path,num2str(i),'.jpg']);
        end
    end
end