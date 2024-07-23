function [node, element] = ...
    quad_mesh(Lx, Ly, nelemX, nelemY, ...
    elementType)

deltaX = Lx/nelemX;
deltaY = Ly/nelemY;
switch elementType
    case 4        
        nodesX = nelemX+1;
        nodesY = nelemY+1;
        node = [];
        for j = 1:nodesY
            for i = 1:nodesX
                x = (i-1)*deltaX; y = (j-1)*deltaY;
                node = [node; x y];
            end
        end
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1 = i+(j-1)*nodesX;
                i2 = i1+1;
                i3 = i2+nodesX;
                i4 = i1+nodesX;
                element = [element; i1 i2 i3 i4];
            end
        end
        
    case 8 
        nodesX = nelemX*2+1;
        nodesY = nelemY*2+1;
        node = [];
        for j = 1:nodesY
            y = (j-1)*deltaY/2;
            if mod(j,2) 
                nodesX = nelemX*2+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX/2;
                    node = [node; x y];
                end
            else 
                nodesX = nelemX+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX;
                    node = [node; x y];
                end
            end
        end
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1 = 1 + 2*(i-1) + (j-1)*round(1.5*nodesX);
                i2 = i1 + 2;
                i3 = i2 + round(1.5*nodesX);
                i4 = i1 + round(1.5*nodesX);
                i5 = i1 + 1;
                i6 = i1 + nodesX+2-i;
                i7 = i4 + 1;
                i8 = i6 - 1;
                element = [element; i1 i2 i3 i4 i5 i6 i7 i8];
            end
        end
        
    case 12 
        nodesX = nelemX*3+1;
        nodesY = nelemY*3+1;
        node = [];
        k = 1;
        for j = 1:nodesY
            y = (j-1)*deltaY/3;
            if j == k
                nodesX = nelemX*3+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX/3;
                    node = [node; x y];
                    k = k + (3/nodesX);
                end
            else 
                nodesX = nelemX+1;
                for i = 1:nodesX
                    x = (i-1)*deltaX;
                    node = [node; x y];
                end
            end
            k = round(k);
        end
        element = [];
        for j = 1:nelemY
            for i = 1:nelemX
                i1  = 1 + 3*(i-1) + (j-1)*(2*nodesX) - 2*(j-1);
                i2  = i1 + 3;
                i3  = i1 + nodesX + (nelemX*2+2) - 2*(i-2.5) + 2*(i-1);
                i4  = i3 - 3;
                i5  = i1 + 1;
                i6  = i5 + 1;
                i7  = i1 + nodesX + 1 - 2*(i-1);
                i8  = i1 + nodesX + (nelemX+2) - 2*(i-1);
                i9  = i3 - 1;
                i10 = i9 - 1;
                i11 = i8 - 1;
                i12 = i7 - 1;
                element = [element; i1 i2 i3 i4 i5 i6 i7 i8 i9 i10 i11 i12];
            end
        end
end
end