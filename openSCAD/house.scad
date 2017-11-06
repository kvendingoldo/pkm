translate([0,0,5]) {
    rotate(a=[0,180,0]) {
        color("red")
        linear_extrude(height = 2, scale=3) 
        square ([3,2],center = true);
    }
};

translate([0,0,0]) {
    color("gray")
    linear_extrude(height = 3)
    square ([9,6],center = true);
};

translate([-0.75,3,0]) {
    cube([1.5,0.5,2.5]);
};


translate([2.5,3,1]) {
    rotate(a=[90,0,0]) {
        polygon(points=[[0,0],[0,1],[0.5,1.5],[1,1],[1,0]]);
    }
}

translate([-3,3,1]) {
    rotate(a=[90,0,0]) {
        polygon(points=[[0,0],[0,1],[0.5,1.5],[1,1],[1,0]]);
    }
}
