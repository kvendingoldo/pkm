echo(version=version());

font = "Liberation Sans";

cube_size = 12;
letter_size = 10;
letter_height = 1;

o = cube_size / 2 - letter_height / 2;

module letter(l) {
    linear_extrude(height = letter_height) {
        text(l, size = letter_size, font = font, halign = "center", valign = "center", $fn = 16);
    }
}

translate([0,-24,0])

difference() {
    union() {
        color("grey") cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 0, 90]) letter("S");
    }
    translate([0, 0, o])  letter("\u263A");
}

translate([0,-12,7])

difference() {
    union() {
        color("grey") rotate([50, 0, 330]) cube(cube_size, center = true);
    }
}


translate([10,-15,-6])
difference () {
    union() {
        translate([o, 0, 0]) rotate([0, 0, 310]) letter("H");
    }
}

translate([0,0,0])

difference() {
    union() {
        color("grey") rotate([90, 0, 0]) cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 0, 90]) letter("A");
    }
}

translate([0,12,8])

difference() {
    union() {
        color("grey") rotate([45, 0, 15]) cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 45, 105]) letter("R");
    }
}

translate([0,24,0])

difference() {
    union() {
        color("grey") rotate([90, 0, 0]) cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 0, 90]) letter("O");
    }
}

translate([0,34,10])

difference() {
    union() {
        color("grey") rotate([90, 0, 0]) cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 0, 90]) letter("V");
    }
}

translate([0,43,0])

difference() {
    union() {
        color("grey") rotate([90, 0, 0]) cube(cube_size, center = true);
        translate([o, 0, 0]) rotate([90, 0, 90]);
    }
}
