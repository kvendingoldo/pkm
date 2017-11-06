module wing(scale=10) {
    rotate([90, 0, 0])
    linear_extrude(height = scale*0.02, center = true)
    polygon(points=[[0, 0], [0, scale], [scale*0.7, scale*0.7], [scale, 0]]);
}

module rocket(scale=10) {
    rotate([0, 0, 0])
    wing(scale);

    rotate([0, 0, 90])
    wing(scale);

    rotate([0, 0, 180])
    wing(scale);

    rotate([0, 0, 270])
    wing(scale);

    color("blue", 0.8)
    cylinder(h = scale, r1=scale/4, r2=scale/3);

    mid_height = scale/2;
    mid_width  = scale/3.5;
    translate([0, 0, scale])
    color("green",0.8)
    cylinder(h = mid_height, r1=scale/3, r2=mid_width);

    long_height = scale+mid_height;
    translate([0, 0, long_height])
    color("blue", 0.8)
    cylinder(h = scale*1.75, r1=mid_width, r2=scale/3.5);

    cone_height = scale*1.75+long_height;
    translate([0, 0, cone_height])
    color("black",0.8)
    cylinder(h = (scale/2), r1=mid_width, r2=0);

}

rocket(29);
