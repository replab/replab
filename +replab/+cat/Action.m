classdef Action
    
    properties (SetAccess = protected)
        groupCat;
    end
    
    methods
        
        function p1 = leftAction(self, g, p)
            error('Not implemented');
        end
        
        function p1 = rightAction(self, p, g)
            error('Not implemented');
        end
        
    end
    
    methods
        
        function verifyLaws(self, randomG, randomP)
            x = randomG();
            y = randomG();
            xy = self.groupCat.compose(x, y);
            p = randomP();
            id = self.groupCat.identity;
            assertEqual(self.leftAction(id, p), p);
            assertEqual(self.rightAction(p, id), p);
            x_y_p = self.leftAction(x, self.leftAction(y, p));
            xy_p = self.leftAction(xy, p);
            assertEqual(x_y_p, xy_p);
            p_x_y = self.rightAction(self.rightAction(p, x), y);
            p_xy = self.rightAction(p, xy);
            assertEqual(p_x_y, p_xy);
            xInv = self.groupCat.inverse(x);
            assertEqual(self.leftAction(x, p), self.rightAction(p, xInv));
        end
        
    end
    
end
