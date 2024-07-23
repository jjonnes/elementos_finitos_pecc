function draw_mesh(nodeCoordinates,elementNodes,elementType,lineType)
%drawingMesh mesh representation
% Input data required
%
% nodeCoordinates : Cartesian coordinates of the nodes
% elementNodes    : Element connectivity
% elementType     : Type of element (e.g. linear, quadratic, etc.)
% lineType        : Type of line for plotting purposes
%%
% resize the vector to suit 2d and 3d problems (beams)
if size(nodeCoordinates,2) == 2
    nodeCoordinates(end,3) = 0;
end
switch elementType
    case 'Q4'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,:),1),nodeCoordinates(elementNodes(i,:),2),'w','FaceColor','none','LineStyle',lineType,'EdgeColor','k')
        end
    case 'Q8'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,1:4),1),nodeCoordinates(elementNodes(i,1:4),2),'w','FaceColor','none','LineStyle',lineType,'EdgeColor','k')
        end
    case 'Q12'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,1:4),1),nodeCoordinates(elementNodes(i,1:4),2),'w','FaceColor','none','LineStyle',lineType,'EdgeColor','k')
        end
    otherwise
        disp('Element type not available')
end
end