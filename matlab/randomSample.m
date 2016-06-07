function [Y] = randomSample(X,h,w,a,b)
    Y=zeros(a*b,size(X,2));
    for i=1:size(X,2)
        im=reshape(X(:,i),h,w);
        x=randi([1,h-a+1],1,1);
        y=randi([1,w-b+1],1,1);
        Y(:,i)=reshape(im(y:y-1+a,x:x-1+b),a*b,1);
    end
end