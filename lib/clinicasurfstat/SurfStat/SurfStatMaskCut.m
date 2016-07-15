function mask = SurfStatMaskCut( surf );

%Mask that excludes the inter-hemisphere cut.
%
% Usage: mask = SurfStatMaskCut( surf );
%
% surf.coord   = 3 x v matrix of surface coordinates, v=#vertices.
%
% mask         = 1 x v vector, 1=inside, 0=outside, v=#vertices.
%
% It looks in -50<y<50 and -20<z<40, and mask vertices where |x|>thresh,
% where thresh = 1.5 x arg max of a histogram of |x|. 

% In this function, we have a bug, cuz here when you use the function, when
% you define the b, and use hist, maybe you will get two max values, so you
% have to change the interval number by mannual, so that is not good, in
% Alex, code, his method is better, cuz, he just check out every vertex
% which is equal to 0, and make them to be a mask!

f=(abs(surf.coord(2,:))<50 ...
    & abs(surf.coord(3,:)-10)<30 ...
    & abs(surf.coord(1,:))<3);
b=(0:0.01:3);
h=hist(abs(surf.coord(1,f)),b);
t=b(find(h==max(h)))*1.5;
mask=~(abs(surf.coord(1,:))<t & f);

return
end
