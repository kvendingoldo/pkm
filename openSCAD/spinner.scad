height = 7;
diameter = 23;
wings = 4;
shell = 4;

union() {
    bearingCylinder(height,diameter);
    translate([2,-8,0]) wings();
    degrees = 360 / wings;
    for (i = [1:wings]) {
        rotate(degrees * i) translate([2,-8,0]) wings();
    }
}

module bearingCylinder (height = height, diameter = diameter) {
    inner = diameter + 1.5;
    outer = diameter + shell;
    difference() {
        //https://doc.cgal.org/latest/Minkowski_sum_3/index.html
        minkowski() {
            $fn = 360;
            cylinder(height-4, r = outer/2);
            $fn = 20;
            sphere(2);
        }
        $fn = 360;
        translate([0,0,-2]) cylinder(height, r = inner/2);
    }
}
module wings (height = height, diameter = diameter) {
    union() {
        difference() {
            outer = diameter + shell;
            minkowski() {
                cube([25,20,3]);
                $fn = 20;
                sphere(2);
            }
            $fn = 360;
            translate([-2,8,-2]) cylinder(height, r = outer/2);
            translate([25,8,-2]) cylinder(height, r = outer/2);
        }
        translate([25,8,0]) bearingCylinder(height, diameter);
    }
}
