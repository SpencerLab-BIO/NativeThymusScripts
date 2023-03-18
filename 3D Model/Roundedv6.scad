//Circular Thymus Holder
//Glass CoverSlip Diameter = 3mm
$fn=1000;

//Height of ring
//height = .5;
height = 0.3;
height2 = 0.5;

//difference = .45

//1.5mm outer edge
//outer edge
//diameter = 1.25;
//diameter2 = 1.25;
//diameter = 1.75
diameter = 1.75;
diameter2 = diameter;

//inner edge
reduced_diameter = diameter - 0.45;

//box height
l_c2 = 6;

//hole
hole=3;


module RingSize(){
    translate([2,-3,0.75]){
        rotate([0,0,90]){
            linear_extrude(0.75)
            text(str(20*reduced_diameter),size=4);
        }
    }
}

//base circle
module BaseRing(){
    translate([0,diameter-1.75,0]){
        difference(){
            cylinder(height2,r1=diameter,r2=diameter2);
            cylinder(height2,r=reduced_diameter);
        };
    }
}

//Extension Arm
module ExtArm(){
    widthEA = 2;
    widthHalf = widthEA / 2;
    translate([-widthHalf,-10,3]){
        rotate([-20,0,0]){
            cube([widthEA,9,height2]);
        };
    };
};

module ExtArm2(){
    widthEA = 2;
    widthHalf = widthEA / 2;
    translate([0,-5.7,1.6]){
        rotate([73,0,0]){
            linear_extrude(height=9,center=true,scale=2)
            square([widthEA,0.5],center=true);
        };
    };
};

//End Piece
module EndPiece(){
    translate([0,0.2,0.4]){
        //Extension arm
        rotate([-10,0,0]){
          translate([0,-17.5,.75]){
            cube([4,15,1.5],center=true);
            RingSize();
          };  
        };

        //connection post
        difference(){
            translate([0,-27,7.4]){
                cube([l_c2,l_c2,l_c2],center=true);
                //RingSize();
            };
            translate([0,-26.5,7.4]){
                   rotate([90,0,0]){
                        cylinder(6,d=hole);
                };
            };
        };
    };
};

//Height Marker
module HeightMarker(){
    translate([0,-6.3,1.8]){
        cube([5,1,.45],center=true);
    }
};

Dist1 = 3.5;
module ObjectiveModel(){
    //diameterO1 = 6.2;
    //diameterO2 = 9.7;
    diameterO1 = 7;
    diameterO2 = 14.6;
    diameterO3 = 15.7;
    
    heightO1 = 41.5-37.7;
    heightO2 = 2.5;
    
    cylinder(h=heightO1,d1=diameterO1,d2=diameterO2);
    translate([0,0,heightO1]){
        cylinder(h=heightO2,d1=diameterO2,d2 = diameterO3);
    };
   
    translate([0,0,-Dist1]){
        cube([0.3,0.3,0.1],center=true);
    };
};

Dist2 = 2;
module ObjectiveModel2(){
    diameterO1 = 6.2;
    diameterO2 = 9.7;
    diameterO3 = 15.7;
    
    heightO1 = 1.3;
    heightO2 = 2.5;
    
    cylinder(h=heightO1,d1=diameterO1,d2=diameterO2);
    translate([0,0,heightO1]){
        cylinder(h=heightO2,d1=diameterO2,d2 = diameterO3);
    };
   
    translate([0,0,-Dist2]){
        cube([0.3,0.3,0.1],center=true);
    };
};

BaseRing();
ExtArm2();
EndPiece();
HeightMarker();

translate([0,0,Dist1]){
    //ObjectiveModel();
};
translate([0,0,Dist2]){
    //ObjectiveModel2();
};

